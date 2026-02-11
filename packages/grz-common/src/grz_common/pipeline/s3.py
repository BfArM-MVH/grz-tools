import hashlib
import io
import logging
from concurrent.futures import ThreadPoolExecutor
from typing import Any

log = logging.getLogger(__name__)


class S3Downloader(io.BufferedIOBase):
    def __init__(self, s3_client: Any, bucket: str, key: str):
        self.response = s3_client.get_object(Bucket=bucket, Key=key)
        self.stream = self.response["Body"]
        self.length = self.response.get("ContentLength", 0)

    def read(self, size: int = -1) -> bytes:
        return self.stream.read(size)

    def close(self) -> None:
        self.stream.close()
        super().close()


class S3MultipartUploader:
    def __init__(self, s3_client: Any, bucket: str, key: str, part_size: int = 64 * 1024 * 1024, max_threads: int = 4):
        self.s3 = s3_client
        self.bucket = bucket
        self.key = key
        self.part_size = part_size
        self.executor = ThreadPoolExecutor(max_workers=max_threads)
        self._upload_id: str | None = None
        self._parts: list[dict[str, Any]] = []
        self._error: BaseException | None = None

    def upload(self, input_stream: io.BufferedIOBase) -> None:
        try:
            response = self.s3.create_multipart_upload(Bucket=self.bucket, Key=self.key)
            self._upload_id = response["UploadId"]

            futures = []
            part_number = 1

            while True:
                if self._error:
                    raise self._error

                chunk = input_stream.read(self.part_size)
                if not chunk:
                    break

                future = self.executor.submit(self._upload_part, self._upload_id, part_number, chunk)
                futures.append(future)
                part_number += 1

                # backpressure
                while len(futures) > 4:
                    self._collect_part(futures.pop(0))

            # collect remaining parts
            for f in futures:
                self._collect_part(f)

            if self._error:
                raise self._error

            self._parts.sort(key=lambda x: x["PartNumber"])
            expected_etag = self._calculate_multipart_etag(self._parts)
            complete_response = self.s3.complete_multipart_upload(
                Bucket=self.bucket,
                Key=self.key,
                UploadId=self._upload_id,
                MultipartUpload={"Parts": [{"PartNumber": p["PartNumber"], "ETag": p["ETag"]} for p in self._parts]},
            )

            server_etag = complete_response.get("ETag", "").strip('"')
            if expected_etag and server_etag != expected_etag:
                raise OSError(f"S3 Upload Corrupt: ETag mismatch! Expected {expected_etag}, Got {server_etag}")

            log.info(f"Upload verified: {self.key} (ETag: {server_etag})")

        except Exception as e:
            self.abort()
            raise e
        finally:
            self.executor.shutdown(wait=False)

    def _upload_part(self, uid: str, part_num: int, data: bytes) -> dict[str, Any]:
        local_md5 = hashlib.md5(data, usedforsecurity=False).digest()

        resp = self.s3.upload_part(Bucket=self.bucket, Key=self.key, UploadId=uid, PartNumber=part_num, Body=data)
        return {
            "PartNumber": part_num,
            "ETag": resp["ETag"],
            "local_md5": local_md5,  # Store raw bytes for final calculation
        }

    def _calculate_multipart_etag(self, sorted_parts: list[dict[str, Any]]) -> str:
        """
        Reconstructs the S3 Multipart ETag: MD5(Concat(PartMD5s)) - NumParts
        """
        md5_digests = [p["local_md5"] for p in sorted_parts if "local_md5" in p]

        if not md5_digests:
            return ""

        combined = b"".join(md5_digests)
        combined_md5 = hashlib.md5(combined, usedforsecurity=False).hexdigest()
        return f"{combined_md5}-{len(md5_digests)}"

    def _collect_part(self, future) -> None:
        try:
            self._parts.append(future.result())
        except Exception as e:
            self._error = e

    def abort(self) -> None:
        if self._upload_id:
            try:
                self.s3.abort_multipart_upload(Bucket=self.bucket, Key=self.key, UploadId=self._upload_id)
            except:
                pass
