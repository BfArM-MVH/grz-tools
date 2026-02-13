import contextlib
import hashlib
import io
import logging
import math
from concurrent.futures import ThreadPoolExecutor
from typing import Any

log = logging.getLogger(__name__)


MULTIPART_DEFAULT_PART_SIZE = 8 * 1024 * 1024  # 8MiB
MULTIPART_MAX_PARTS = 1000  # Ceph S3 limit
MULTIPART_MIN_PART_SIZE = 5 * 1024 * 1024  # S3 Standard Min (5MB)


def calculate_s3_part_size(file_size: int | None, user_part_size: int | None = None) -> int:
    """
    Calculate the appropriate part size for a multipart upload.
    Ensures the number of parts stays within limits (Ceph/Bluestore/Quobyte: 1000).
    """
    if user_part_size is not None:
        return user_part_size

    if file_size is None or file_size <= 0:
        return MULTIPART_DEFAULT_PART_SIZE

    # If default part size would result in too many parts, increase it
    if file_size / MULTIPART_DEFAULT_PART_SIZE > MULTIPART_MAX_PARTS:
        optimal = math.ceil(file_size / MULTIPART_MAX_PARTS)
        return max(optimal, MULTIPART_MIN_PART_SIZE)

    return max(MULTIPART_DEFAULT_PART_SIZE, MULTIPART_MIN_PART_SIZE)


class S3Downloader(io.BufferedIOBase):
    def __init__(self, s3_client: Any, bucket: str, key: str):
        self.response = s3_client.get_object(Bucket=bucket, Key=key)
        self.stream = self.response["Body"]
        self.length = self.response.get("ContentLength", 0)

    def read(self, size: int | None = -1) -> bytes:
        try:
            if size == -1 or size is None:
                return self.stream.read()
            return self.stream.read(size)
        except Exception as e:
            raise OSError(f"S3 Read Error: {e}") from e

    def close(self) -> None:
        if hasattr(self, "stream"):
            self.stream.close()
        super().close()


class S3MultipartUploader(io.BufferedIOBase):
    """
    A write-only stream that uploads data to S3 using Multipart Upload.
    """

    def __init__(
        self,
        s3_client: Any,
        bucket: str,
        key: str,
        part_size: int | None = None,
        max_threads: int = 4,
    ):
        super().__init__()
        self.s3 = s3_client
        self.bucket = bucket
        self.key = key
        self.part_size = calculate_s3_part_size(None, part_size)
        self.max_threads = max_threads

        self._executor: ThreadPoolExecutor | None = None
        self._upload_id: str | None = None
        self._parts: list[dict[str, Any]] = []
        self._futures: list[Any] = []
        self._buffer = bytearray()
        self._part_number = 1
        self._closed = False

    def writable(self) -> bool:
        return True

    def __enter__(self):
        self._start_multipart_upload()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_type:
            self.abort()
        else:
            self.close()

    def _start_multipart_upload(self):
        if self._upload_id:
            return
        log.info(f"S3Uploader: Starting upload to s3://{self.bucket}/{self.key}")
        try:
            resp = self.s3.create_multipart_upload(Bucket=self.bucket, Key=self.key)
            self._upload_id = resp["UploadId"]
            self._executor = ThreadPoolExecutor(max_workers=self.max_threads)
        except Exception:
            self._cleanup()
            raise

    def write(self, b: bytes) -> int:
        """
        Write bytes to the buffer. Uploads parts when buffer is full.
        """
        if self._closed:
            raise ValueError("I/O operation on closed file.")

        if not self._upload_id:
            self._start_multipart_upload()

        self._check_futures()

        self._buffer.extend(b)
        while len(self._buffer) >= self.part_size:
            chunk = self._buffer[: self.part_size]
            del self._buffer[: self.part_size]
            self._submit_part(chunk, self._part_number)
            self._part_number += 1

        return len(b)

    def close(self) -> None:
        """
        Flush remaining buffer, wait for threads, and complete upload.
        """
        if self._closed:
            return

        self._closed = True

        try:
            # upload remaining data
            if self._buffer:
                self._submit_part(bytes(self._buffer), self._part_number)
                self._buffer.clear()

            if self._futures:
                for f in self._futures:
                    self._parts.append(f.result())

            if self._upload_id:
                self._parts.sort(key=lambda x: x["PartNumber"])
                self._complete_upload()

        except Exception as e:
            log.error(f"Upload failed: {e}")
            self.abort()
            raise e
        finally:
            self._cleanup()

    def abort(self) -> None:
        """Abort the multipart upload on S3."""
        if self._upload_id:
            with contextlib.suppress(Exception):
                self.s3.abort_multipart_upload(Bucket=self.bucket, Key=self.key, UploadId=self._upload_id)
        self._cleanup()

    def _submit_part(self, data: bytes, part_num: int):
        if not self._executor:
            raise RuntimeError("Executor not initialized")

        future = self._executor.submit(self._upload_part, self._upload_id, part_num, data)
        self._futures.append(future)

    def _upload_part(self, uid: str, part_num: int, data: bytes) -> dict[str, Any]:
        local_md5 = hashlib.md5(data, usedforsecurity=False).digest()
        resp = self.s3.upload_part(
            Bucket=self.bucket,
            Key=self.key,
            UploadId=uid,
            PartNumber=part_num,
            Body=data,
        )
        return {"PartNumber": part_num, "ETag": resp["ETag"], "local_md5": local_md5}

    def _complete_upload(self):
        expected = self._calc_etag(self._parts)
        parts_payload = [{"PartNumber": p["PartNumber"], "ETag": p["ETag"]} for p in self._parts]

        if not parts_payload:
            empty_part = self._upload_part(self._upload_id, 1, b"")
            parts_payload.append({"PartNumber": 1, "ETag": empty_part["ETag"]})

        complete = self.s3.complete_multipart_upload(
            Bucket=self.bucket,
            Key=self.key,
            UploadId=self._upload_id,
            MultipartUpload={"Parts": parts_payload},
        )

        server_etag = complete.get("ETag", "").strip('"')
        if expected and server_etag != expected:
            raise OSError(f"ETag mismatch! Exp: {expected}, Got: {server_etag}")

    def _calc_etag(self, parts: list[dict[str, Any]]) -> str:
        digests = [p["local_md5"] for p in parts if "local_md5" in p]
        if not digests:
            return ""
        combined = b"".join(digests)
        combined_hash = hashlib.md5(combined, usedforsecurity=False).hexdigest()
        return f"{combined_hash}-{len(digests)}"

    def _check_futures(self):
        """Check if any background tasks failed."""
        self._futures = [f for f in self._futures if not f.done() or f.result()]

    def _cleanup(self):
        if self._executor:
            self._executor.shutdown(wait=False)
            self._executor = None
