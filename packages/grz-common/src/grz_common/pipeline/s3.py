"""S3 source and sink pipeline stages for streaming access."""

from __future__ import annotations

import contextlib
import hashlib
import io
import logging
import math
from concurrent.futures import Future, ThreadPoolExecutor
from typing import Any

import botocore.exceptions

from .base import PipelineError, StreamSink, StreamSource
from .constants import (
    MULTIPART_DEFAULT_PART_SIZE,
    MULTIPART_MAX_PARTS,
    MULTIPART_MIN_PART_SIZE,
)

log = logging.getLogger(__name__)


class S3Downloader(StreamSource):
    """Streaming source that downloads from S3."""

    def __init__(
        self,
        s3_client: Any,
        bucket: str,
        key: str,
        name: str | None = None,
    ):
        """
        Initialize the S3 downloader.

        :param s3_client: Boto3 S3 client
        :param bucket: S3 bucket name
        :param key: S3 object key
        :param name: Stage name for logging
        """
        super().__init__(name or "S3Downloader")
        self._s3_client = s3_client
        self._bucket = bucket
        self._key = key
        self._body: Any = None
        self._content_length: int | None = None
        self._etag: str | None = None
        self._bytes_read = 0

    def initialize(self, context: Any) -> None:
        """Open the S3 object for streaming."""
        super().initialize(context)

        try:
            # get object metadata
            head = self._s3_client.head_object(Bucket=self._bucket, Key=self._key)
            self._content_length = head["ContentLength"]
            self._etag = head.get("ETag", "").strip('"')

            # open streaming body
            response = self._s3_client.get_object(Bucket=self._bucket, Key=self._key)
            self._body = response["Body"]

            self._log.debug(f"Opened S3 object s3://{self._bucket}/{self._key} ({self._content_length} bytes)")

        except botocore.exceptions.ClientError as e:
            error_code = e.response.get("Error", {}).get("Code", "")
            if error_code == "404":
                raise PipelineError(
                    f"S3 object not found: s3://{self._bucket}/{self._key}",
                    self.name,
                    e,
                ) from e
            raise PipelineError(
                f"Failed to open S3 object: {e}",
                self.name,
                e,
            ) from e

    def read(self, size: int = -1) -> bytes:
        """
        Read data from the S3 object.

        :param size: Number of bytes to read (-1 for all)
        :returns: Bytes read, empty when exhausted
        """
        if self._body is None:
            return b""

        data = self._body.read() if size < 0 else self._body.read(size)

        if data:
            self._bytes_read += len(data)

        return data

    def _close_body(self) -> None:
        """Close the S3 stream body if open."""
        if self._body is not None:
            with contextlib.suppress(Exception):
                self._body.close()
            self._body = None

    def finalize(self) -> None:
        """Close the S3 stream and record stats."""
        self._close_body()
        self.context.bytes_read = self._bytes_read
        self.context.metadata["s3_source_key"] = self._key
        self.context.metadata["s3_source_etag"] = self._etag

    def abort(self) -> None:
        """Close the S3 stream on abort."""
        self._close_body()

    @property
    def content_length(self) -> int | None:
        """Return the content length of the S3 object."""
        return self._content_length

    @property
    def etag(self) -> str | None:
        """Return the ETag of the S3 object."""
        return self._etag

    @property
    def bytes_read(self) -> int:
        """Return total bytes read."""
        return self._bytes_read


def calculate_part_size(
    file_size: int,
    default_part_size: int = MULTIPART_DEFAULT_PART_SIZE,
    max_parts: int = MULTIPART_MAX_PARTS,
    min_part_size: int = MULTIPART_MIN_PART_SIZE,
) -> int:
    """
    Calculate the appropriate part size for a multipart upload.

    Ensures the number of parts stays within limits (Ceph/Bluestore/Quobyte: 1000).

    :param file_size: Size of the file to upload in bytes
    :param default_part_size: Default part size if file is small enough
    :param max_parts: Maximum number of parts allowed
    :param min_part_size: Minimum part size (S3 requires 5MB)
    :returns: Part size to use in bytes
    """
    if file_size <= 0:
        return default_part_size

    # if default part size would result in too many parts, increase it
    if file_size / default_part_size > max_parts:
        calculated_size = math.ceil(file_size / max_parts)
        return max(calculated_size, min_part_size)

    return max(default_part_size, min_part_size)


class S3MultipartUploader(StreamSink):
    """
    Streaming sink that uploads to S3 using multipart upload.

    Computes MD5 checksums for each part and verifies the final ETag
    to ensure data integrity during upload.
    """

    def __init__(  # noqa: PLR0913
        self,
        s3_client: Any,
        bucket: str,
        key: str,
        part_size: int = MULTIPART_DEFAULT_PART_SIZE,
        max_concurrent_uploads: int = 1,
        expected_size: int | None = None,
        name: str | None = None,
    ):
        """
        Initialize the S3 multipart uploader.

        :param s3_client: Boto3 S3 client
        :param bucket: S3 bucket name
        :param key: S3 object key
        :param part_size: Size of each part (minimum 5MB, default 256MB)
        :param max_concurrent_uploads: Maximum concurrent part uploads (1 = sequential)
        :param expected_size: Expected total size for dynamic part size calculation
        :param name: Stage name for logging
        """
        super().__init__(name or "S3MultipartUploader")
        self._s3_client = s3_client
        self._bucket = bucket
        self._key = key
        self._max_concurrent = max(1, max_concurrent_uploads)

        # calculate part size based on expected file size to stay within part limits
        if expected_size is not None:
            self._part_size = calculate_part_size(expected_size, part_size)
        else:
            self._part_size = max(part_size, MULTIPART_MIN_PART_SIZE)

        self._upload_id: str | None = None
        self._parts: list[dict[str, Any]] = []
        self._part_md5s: list[bytes] = []  # local MD5 digests for ETag verification
        self._pending_futures: list[Future[dict[str, Any]]] = []
        self._buffer = io.BytesIO()
        self._current_part_number = 1
        self._total_bytes = 0
        self._completed = False
        self._executor: ThreadPoolExecutor | None = None
        self._upload_error: Exception | None = None

    def initialize(self, context: Any) -> None:
        """Initialize the multipart upload."""
        super().initialize(context)

        # create thread pool for concurrent uploads
        if self._max_concurrent > 1:
            self._executor = ThreadPoolExecutor(max_workers=self._max_concurrent)

        try:
            response = self._s3_client.create_multipart_upload(
                Bucket=self._bucket,
                Key=self._key,
            )
            self._upload_id = response["UploadId"]
            self._log.debug(
                f"Initiated multipart upload to s3://{self._bucket}/{self._key} "
                f"(upload_id: {self._upload_id}, concurrent: {self._max_concurrent})"
            )
        except Exception as e:
            self._shutdown_executor()
            raise PipelineError(
                f"Failed to initiate multipart upload: {e}",
                self.name,
                e,
            ) from e

    def _shutdown_executor(self) -> None:
        """Shutdown the thread pool executor."""
        if self._executor is not None:
            self._executor.shutdown(wait=False, cancel_futures=True)
            self._executor = None

    def write(self, data: bytes) -> int:
        """
        Write data to the upload buffer.

        Automatically uploads parts when buffer reaches part_size.

        :param data: Bytes to write
        :returns: Number of bytes written
        :raises PipelineError: If a concurrent upload has failed
        """
        if not data:
            return 0

        # check if any concurrent upload has failed
        if self._upload_error is not None:
            raise PipelineError(
                f"Upload failed: {self._upload_error}",
                self.name,
                self._upload_error,
            )

        # collect completed futures to manage memory and detect errors early
        self._collect_completed_futures()

        self._buffer.write(data)
        self._total_bytes += len(data)

        # upload complete parts
        while self._buffer.tell() >= self._part_size:
            self._upload_part()

        return len(data)

    def _collect_completed_futures(self) -> None:
        """Collect results from completed upload futures."""
        still_pending: list[Future[dict[str, Any]]] = []
        for future in self._pending_futures:
            if future.done():
                try:
                    part_info = future.result()
                    self._parts.append(part_info)
                except Exception as e:
                    self._upload_error = e
                    self._log.error(f"Part upload failed: {e}")
            else:
                still_pending.append(future)
        self._pending_futures = still_pending

    def _upload_part(self) -> None:
        """Upload the current buffer as a part (sync or async via thread pool)."""
        if self._buffer.tell() == 0:
            return

        if self._upload_id is None:
            raise RuntimeError("Upload not initialized")

        self._buffer.seek(0)
        part_data = self._buffer.read(self._part_size)
        remaining = self._buffer.read()

        part_number = self._current_part_number
        self._current_part_number += 1

        if self._executor is not None:
            # submit to thread pool for concurrent upload
            future = self._executor.submit(self._upload_part_sync, part_data, part_number)
            self._pending_futures.append(future)
            self._log.debug(f"Submitted part {part_number} for upload ({len(part_data)} bytes)")
        else:
            # synchronous upload
            try:
                part_info = self._upload_part_sync(part_data, part_number)
                self._parts.append(part_info)
            except Exception as e:
                raise PipelineError(
                    f"Failed to upload part {part_number}: {e}",
                    self.name,
                    e,
                ) from e

        # keep remaining data in buffer
        self._buffer = io.BytesIO()
        if remaining:
            self._buffer.write(remaining)

    def _upload_part_sync(self, part_data: bytes, part_number: int) -> dict[str, Any]:
        """
        Synchronously upload a single part to S3.

        Computes MD5 checksum locally and verifies against the returned ETag.

        :param part_data: The part data to upload
        :param part_number: The part number (1-indexed)
        :returns: Dict with PartNumber, ETag, and local_md5
        :raises PipelineError: If upload fails or checksum mismatch
        """
        # compute MD5 checksum locally
        local_md5 = hashlib.md5(part_data, usedforsecurity=False)
        local_md5_hex = local_md5.hexdigest()

        response = self._s3_client.upload_part(
            Bucket=self._bucket,
            Key=self._key,
            UploadId=self._upload_id,
            PartNumber=part_number,
            Body=part_data,
        )

        # S3 ETag is the MD5 hex wrapped in quotes
        s3_etag = response["ETag"].strip('"')

        # verify checksum matches
        if s3_etag != local_md5_hex:
            raise PipelineError(
                f"Part {part_number} checksum mismatch: local={local_md5_hex}, S3={s3_etag}",
                self.name,
            )

        self._log.debug(f"Uploaded part {part_number} ({len(part_data)} bytes, ETag: {s3_etag}, verified)")

        return {
            "PartNumber": part_number,
            "ETag": response["ETag"],
            "local_md5": local_md5.digest(),
        }

    def finalize(self) -> None:
        """Complete the multipart upload and verify final ETag."""
        if self._completed:
            return

        if self._upload_id is None:
            raise RuntimeError("Upload not initialized")

        try:
            # upload any remaining data in buffer as final part
            if self._buffer.tell() > 0:
                self._buffer.seek(0)
                remaining_data = self._buffer.read()
                part_number = self._current_part_number

                if self._executor is not None:
                    # submit final part to thread pool
                    future = self._executor.submit(self._upload_part_sync, remaining_data, part_number)
                    self._pending_futures.append(future)
                else:
                    # synchronous upload
                    part_info = self._upload_part_sync(remaining_data, part_number)
                    self._parts.append(part_info)

            # wait for all pending uploads to complete
            self._wait_for_pending_uploads()

            # check for upload errors
            if self._upload_error is not None:
                raise PipelineError(
                    f"One or more part uploads failed: {self._upload_error}",
                    self.name,
                    self._upload_error,
                )

            # complete the upload and verify ETag
            if self._parts:
                # sort parts by part number (required by S3)
                sorted_parts = sorted(self._parts, key=lambda p: p["PartNumber"])

                # compute expected ETag from local MD5s
                expected_etag = self._compute_multipart_etag(sorted_parts)

                # prepare parts for S3 (only PartNumber and ETag)
                s3_parts = [{"PartNumber": p["PartNumber"], "ETag": p["ETag"]} for p in sorted_parts]

                response = self._s3_client.complete_multipart_upload(
                    Bucket=self._bucket,
                    Key=self._key,
                    UploadId=self._upload_id,
                    MultipartUpload={"Parts": s3_parts},
                )

                # verify final ETag - critical for data integrity
                actual_etag = response.get("ETag", "").strip('"')
                if expected_etag and actual_etag and expected_etag != actual_etag:
                    raise PipelineError(
                        f"Final ETag mismatch: expected={expected_etag}, actual={actual_etag}. "
                        "Upload data may be corrupted.",
                        self.name,
                    )

                self._log.debug(
                    f"Completed multipart upload to s3://{self._bucket}/{self._key} "
                    f"({len(self._parts)} parts, {self._total_bytes} bytes, ETag: {actual_etag})"
                )

                self.context.metadata["s3_etag"] = actual_etag
            else:
                # no parts uploaded: abort and create empty object
                self._s3_client.abort_multipart_upload(
                    Bucket=self._bucket,
                    Key=self._key,
                    UploadId=self._upload_id,
                )
                self._s3_client.put_object(
                    Bucket=self._bucket,
                    Key=self._key,
                    Body=b"",
                )
                self._log.debug(f"Created empty object s3://{self._bucket}/{self._key}")

            self._completed = True
            self.context.bytes_written = self._total_bytes
            self.context.metadata["s3_target_key"] = self._key
            self.context.metadata["s3_parts_uploaded"] = len(self._parts)

        except PipelineError:
            self.abort()
            raise
        except Exception as e:
            self.abort()
            raise PipelineError(
                f"Failed to complete multipart upload: {e}",
                self.name,
                e,
            ) from e
        finally:
            self._shutdown_executor()

    def _wait_for_pending_uploads(self) -> None:
        """Wait for all pending upload futures to complete and collect results."""
        for future in self._pending_futures:
            try:
                part_info = future.result()
                self._parts.append(part_info)
            except Exception as e:
                if self._upload_error is None:
                    self._upload_error = e
                self._log.error(f"Part upload failed: {e}")
        self._pending_futures = []

    def _compute_multipart_etag(self, sorted_parts: list[dict[str, Any]]) -> str:
        """
        Compute the expected multipart ETag from local MD5 digests.

        S3 multipart ETag format: MD5(concat(part_md5s))-num_parts
        This applies even for single-part multipart uploads.

        :param sorted_parts: Parts sorted by PartNumber, each with 'local_md5' key
        :returns: Expected ETag string
        """
        # collect MD5 digests from parts
        md5_digests = []
        for part in sorted_parts:
            if "local_md5" in part:
                md5_digests.append(part["local_md5"])

        if not md5_digests:
            return ""

        # multipart ETag is always: MD5(concat(part_md5s))-num_parts
        # even for single-part uploads done via multipart API
        combined = b"".join(md5_digests)
        combined_md5 = hashlib.md5(combined, usedforsecurity=False).hexdigest()
        return f"{combined_md5}-{len(md5_digests)}"

    def abort(self) -> None:
        """Abort the multipart upload and cancel pending uploads."""
        if self._completed or self._upload_id is None:
            return

        # cancel pending uploads
        self._shutdown_executor()

        try:
            self._s3_client.abort_multipart_upload(
                Bucket=self._bucket,
                Key=self._key,
                UploadId=self._upload_id,
            )
            self._log.debug(f"Aborted multipart upload {self._upload_id}")
        except Exception as e:
            self._log.warning(f"Failed to abort multipart upload: {e}")

        self._completed = True

    @property
    def bytes_written(self) -> int:
        """Return total bytes written."""
        return self._total_bytes

    @property
    def part_count(self) -> int:
        """Return number of parts uploaded."""
        return len(self._parts)
