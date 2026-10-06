import base64
import hashlib
import logging
import math
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from contextlib import AbstractContextManager
from typing import Any

from botocore.exceptions import (
    ClientError,
    IncompleteReadError,
    ReadTimeoutError,
    ResponseStreamingError,
)
from grz_common.constants import MULTIPART_DEFAULT_PART_SIZE, MULTIPART_MAX_PARTS, MULTIPART_MIN_PART_SIZE
from grz_common.exceptions import (
    DownloadError,
    GrzError,
    MissingObjectError,
    NetworkError,
    UploadError,
)
from grz_common.transfer import s3_errors

from . import Observer, ReadStream, UploadIntegrityError

log = logging.getLogger(__name__)


def calculate_s3_part_size(file_size: int | None, preferred_part_size: int = MULTIPART_DEFAULT_PART_SIZE) -> int:
    """
    Calculate the part size for a multipart upload of ``file_size`` bytes.

    Returns ``preferred_part_size``, at least ``MULTIPART_MIN_PART_SIZE``. If the file would then
    need more than ``MULTIPART_MAX_PARTS`` parts (Ceph/Bluestore/Quobyte: 1000), the part size is
    raised so that it fits.
    """
    part_size = max(preferred_part_size, MULTIPART_MIN_PART_SIZE)
    if file_size is not None and file_size > part_size * MULTIPART_MAX_PARTS:
        part_size = math.ceil(file_size / MULTIPART_MAX_PARTS)
    return part_size


_RESUMABLE_READ_ERRORS = (ReadTimeoutError, IncompleteReadError, ResponseStreamingError, OSError)
"""What reading the body of an S3 object raises when the connection breaks."""


class S3Downloader(ReadStream):
    """Reading from S3 is the Source of the pipeline.

    If the connection breaks while reading, the rest of the object is requested again, from the
    offset reached. The read fails after ``max_resumes`` breaks in a row without data in between.
    ``IfMatch`` makes sure that the rest comes from the same object.
    """

    def __init__(
        self,
        s3_client: Any,
        bucket: str,
        key: str,
        missing_error: type[GrzError] = MissingObjectError,
        max_resumes: int = 5,
    ):
        """
        :param s3_client: boto3 S3 client.
        :param bucket: Name of the bucket.
        :param key: Key of the object.
        :param missing_error: The class for a missing object.
        :param max_resumes: How often the read may break in a row without data in between.
        :raises MissingObjectError: As ``missing_error``, if the object does not exist.
        :raises ConfigurationError: If only a faulty setup causes the error, see :func:`~grz_common.transfer.s3_errors`.
        :raises DownloadError: For any other error of the S3 client.
        """
        # Base first: a stage that raises before it runs is still finalized, and finalizing closes.
        super().__init__()
        self._s3_client = s3_client
        self._bucket = bucket
        self._key = key
        self._missing_error = missing_error
        self._max_resumes = max_resumes
        self._resumes = 0
        self._offset = 0
        # get_object reads only the headers. The body stays on the connection as a StreamingBody,
        # and read() takes it from there in chunks, so the whole object is never in memory.
        self.response = self._get_object()
        self._source = self.response["Body"]
        self.length: int = self.response.get("ContentLength", 0)

    def _get_object(self, **kwargs: Any) -> dict[str, Any]:
        with s3_errors(f"Reading s3://{self._bucket}/{self._key}", DownloadError):
            try:
                return self._s3_client.get_object(Bucket=self._bucket, Key=self._key, **kwargs)
            except ClientError as e:
                # a GET answer carries the real error code, unlike a HEAD answer
                if e.response["Error"]["Code"] == "NoSuchKey":
                    raise self._missing_error(f"s3://{self._bucket}/{self._key} does not exist") from e
                raise

    def read(self, size: int | None = -1) -> bytes:
        while True:
            try:
                chunk = super().read(size)
            except _RESUMABLE_READ_ERRORS as e:
                if self._resumes >= self._max_resumes:
                    raise NetworkError(f"S3 read error: {e}") from e
                self._resumes += 1
                log.warning(f"Reading s3://{self._bucket}/{self._key} broke at byte {self._offset}, resuming: {e}")
                response = self._get_object(Range=f"bytes={self._offset}-", IfMatch=self.response["ETag"])
                self._source = response["Body"]
                status, content_range = response["ResponseMetadata"]["HTTPStatusCode"], response.get("ContentRange", "")
                if status != 206 or not content_range.startswith(f"bytes {self._offset}-"):
                    # a server that ignores Range sends the whole object, which would be appended to the bytes read
                    raise DownloadError(
                        f"Resuming s3://{self._bucket}/{self._key} at byte {self._offset} got status {status} "
                        f"and range {content_range!r} instead of the rest of the object"
                    ) from e
                continue
            except Exception as e:
                raise NetworkError(f"S3 read error: {e}") from e
            if chunk:
                # a long download may break more often than max_resumes, as long as it progresses
                self._resumes = 0
            self._offset += len(chunk)
            return chunk


class S3MultipartUploader(Observer):
    """
    Writing to S3 is a Sink (Observer).
    It buffers data and uploads parts. An empty stream is sent with an empty PUT on close,
    because a multipart upload needs at least one part.
    """

    def __init__(  # noqa: PLR0913, PLR0917
        self,
        s3_client: Any,
        bucket: str,
        key: str,
        part_size: int | None = None,
        max_threads: int = 4,
        content_type: str | None = None,
    ):
        super().__init__()
        self.s3 = s3_client
        self.bucket = bucket
        self.key = key
        self.part_size = part_size or MULTIPART_DEFAULT_PART_SIZE
        self.max_threads = max_threads
        self.content_type = content_type

        self._executor: ThreadPoolExecutor | None = None
        self._upload_id: str | None = None
        self._parts: list[dict[str, Any]] = []
        self._futures: list[Any] = []
        self._buffer = bytearray()
        self._part_number = 1
        self._closed = False

    def __exit__(self, exc_type, exc, tb) -> None:
        # Finish the upload on a clean exit; abort and discard staged parts on error.
        # Returns None (never True), so the original exception is never suppressed.
        if exc is None:
            self.close()
        else:
            self.abort()

    def _start_multipart_upload(self):
        if self._upload_id:
            return
        log.info(f"S3Uploader: Starting upload to s3://{self.bucket}/{self.key}")
        try:
            kwargs: dict[str, Any] = {"Bucket": self.bucket, "Key": self.key}
            if self.content_type:
                kwargs["ContentType"] = self.content_type
            resp = self.s3.create_multipart_upload(**kwargs)
            self._upload_id = resp["UploadId"]
            self._executor = ThreadPoolExecutor(max_workers=self.max_threads)
        except Exception:
            self._cleanup()
            raise

    def _upload_errors(self) -> AbstractContextManager[None]:
        """Raise S3 client errors as an :class:`UploadError`, or as a :class:`ConfigurationError`."""
        return s3_errors(f"Upload to s3://{self.bucket}/{self.key}", UploadError)

    def observe(self, chunk: bytes) -> None:
        """
        Buffers incoming bytes and submits uploads when part_size is reached.
        """
        if self._closed:
            raise ValueError("I/O operation on closed file.")
        if not chunk:
            return

        with self._upload_errors():
            if not self._upload_id:
                self._start_multipart_upload()

            self._check_futures()

            self._buffer.extend(chunk)
            while len(self._buffer) >= self.part_size:
                self._throttle_uploads()
                # the upload takes the buffer itself, and only the overhang is copied into a new one
                part, self._buffer = self._buffer, self._buffer[self.part_size :]
                del part[self.part_size :]
                self._submit_part(part, self._part_number)
                self._part_number += 1

    def _throttle_uploads(self):
        """
        Don't read from disk if the executor queue is full.
        """
        self._check_futures()

        # if we have reached our max concurrency limit, wait for one to finish
        if len(self._futures) >= self.max_threads:
            done, not_done = wait(self._futures, return_when=FIRST_COMPLETED)
            for f in done:
                self._parts.append(f.result())
            self._futures = list(not_done)

    def close(self) -> None:
        """
        Flush remaining buffer, wait for threads, and complete upload.
        """
        if self._closed:
            return

        self._closed = True

        try:
            with self._upload_errors():
                if not self._upload_id:
                    # nothing was written: an empty object needs a PUT
                    self._put_object(self._buffer)
                    self._buffer = bytearray()
                else:
                    # upload remaining data
                    if self._buffer:
                        self._submit_part(self._buffer, self._part_number)
                        self._buffer = bytearray()

                    for f in self._futures:
                        self._parts.append(f.result())

                    self._parts.sort(key=lambda x: x["PartNumber"])
                    self._complete_upload()

        except BaseException as e:
            log.error(f"Upload failed: {e}")
            self.abort()
            raise e
        finally:
            self._cleanup()
            super().close()

    def abort(self) -> None:
        """Abort the multipart upload on S3."""
        if self._executor:
            # a part still in flight could land after the abort and stay in the bucket
            self._executor.shutdown(wait=True, cancel_futures=True)
        if self._upload_id:
            try:
                self.s3.abort_multipart_upload(Bucket=self.bucket, Key=self.key, UploadId=self._upload_id)
            except Exception as e:
                # Don't let a failed abort hide the original error; just log it and move on.
                # A likely cause is missing abort permission, in which case the staged parts
                # stay in the bucket until a lifecycle rule or an admin removes them.
                log.warning(f"Could not abort multipart upload for {self.key}: {e}. Incomplete parts may remain.")
        self._cleanup()
        self._closed = True  # prevent a later close()/finalizer from re-running on an aborted upload

    def _put_object(self, data: bytes | bytearray) -> None:
        hasher = hashlib.md5(data, usedforsecurity=False)
        local_md5_hex = hasher.hexdigest()

        kwargs: dict[str, Any] = {
            "Bucket": self.bucket,
            "Key": self.key,
            "Body": data,
            "ContentMD5": base64.b64encode(hasher.digest()).decode("utf-8"),
        }
        if self.content_type:
            kwargs["ContentType"] = self.content_type
        resp = self.s3.put_object(**kwargs)

        server_etag = resp["ETag"].strip('"')
        if server_etag != local_md5_hex:
            self._delete_mismatched_object()
            raise UploadIntegrityError(
                f"Local checksum does not match remote one! Expected: {local_md5_hex}, Got: {server_etag}",
                stage=self.__class__.__name__,
            )

    def _submit_part(self, data: bytes | bytearray, part_num: int):
        if not self._executor or not self._upload_id:
            raise RuntimeError("Multipart upload not started")

        future = self._executor.submit(self._upload_part, self._upload_id, part_num, data)
        self._futures.append(future)

    def _upload_part(self, uid: str, part_num: int, data: bytes | bytearray) -> dict[str, Any]:
        hasher = hashlib.md5(data, usedforsecurity=False)
        local_md5_bytes = hasher.digest()
        local_md5_hex = hasher.hexdigest()
        b64_md5 = base64.b64encode(local_md5_bytes).decode("utf-8")

        resp = self.s3.upload_part(
            Bucket=self.bucket,
            Key=self.key,
            UploadId=uid,
            PartNumber=part_num,
            Body=data,
            ContentMD5=b64_md5,
        )

        server_etag = resp["ETag"].strip('"')
        if server_etag != local_md5_hex:
            raise UploadIntegrityError(
                f"Local checksum for {part_num} does not match remote one! Expected: {local_md5_hex}, Got: {server_etag}",
                stage=self.__class__.__name__,
            )

        return {"PartNumber": part_num, "ETag": resp["ETag"], "local_md5": local_md5_bytes}

    def _complete_upload(self):
        expected = self._calc_etag(self._parts)
        parts_payload = [{"PartNumber": p["PartNumber"], "ETag": p["ETag"]} for p in self._parts]

        complete = self.s3.complete_multipart_upload(
            Bucket=self.bucket,
            Key=self.key,
            UploadId=self._upload_id,
            MultipartUpload={"Parts": parts_payload},
        )
        # the object exists from here on, so there is no upload left to abort
        self._upload_id = None

        server_etag = complete.get("ETag", "").strip('"')
        if expected and server_etag != expected:
            self._delete_mismatched_object()
            raise UploadIntegrityError(
                f"Final ETag mismatch! Exp: {expected}, Got: {server_etag}", stage=self.__class__.__name__
            )

    def _delete_mismatched_object(self) -> None:
        """Remove the stored object whose content does not match what was sent."""
        try:
            self.s3.delete_object(Bucket=self.bucket, Key=self.key)
        except Exception as e:
            # Don't let a failed delete hide the integrity error; just log it and move on.
            log.warning(f"Could not delete {self.key} after its integrity check failed: {e}. It may still be there.")

    def _calc_etag(self, parts: list[dict[str, Any]]) -> str:
        digests = [p["local_md5"] for p in parts if "local_md5" in p]
        if not digests:
            return ""
        combined = b"".join(digests)
        combined_hash = hashlib.md5(combined, usedforsecurity=False).hexdigest()
        return f"{combined_hash}-{len(digests)}"

    def _check_futures(self):
        """Check if any background tasks failed and collect completed parts."""
        active = []
        for f in self._futures:
            if f.done():
                # .result() raises the exception if the future failed,
                # otherwise it returns the part dictionary we need to keep.
                self._parts.append(f.result())
            else:
                active.append(f)
        self._futures = active

    def _cleanup(self):
        if self._executor:
            self._executor.shutdown(wait=False)
            self._executor = None
