"""Base classes and implementations for pipeline components."""

import contextlib
import hashlib
import io
import logging
import math
import os
import tempfile
import threading
import time
import zlib
from abc import ABCMeta, abstractmethod
from concurrent.futures import ThreadPoolExecutor
from typing import Any

import crypt4gh.header
import crypt4gh.lib
import pysam
from grz_pydantic_models.submission.metadata.v1 import FileType
from nacl.bindings import crypto_aead_chacha20poly1305_ietf_encrypt
from nacl.public import PrivateKey

log = logging.getLogger(__name__)

MULTIPART_DEFAULT_PART_SIZE = 8 * 1024 * 1024  # 8MiB
MULTIPART_MAX_PARTS = 1000  # Ceph S3 limit
MULTIPART_MIN_PART_SIZE = 5 * 1024 * 1024  # S3 Standard Min (5MB)
READ_CHUNK_SIZE = 1 * 1024 * 1024


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


class PipelineError(Exception):
    """Base exception for pipeline errors."""

    def __init__(self, message: str, stage: str | None = None, cause: Exception | None = None):
        self.stage = stage
        self.cause = cause
        super().__init__(f"[{stage}] {message}" if stage else message)


class StreamWrapper(io.BufferedIOBase):
    def __init__(self, source: io.BufferedIOBase):
        self.source = source

    def readable(self) -> bool:
        return True

    def close(self):
        if not self.closed:
            with contextlib.suppress(Exception):
                self.source.close()
            super().close()


class ObserverStream(StreamWrapper, metaclass=ABCMeta):
    @abstractmethod
    def observe(self, chunk: bytes) -> None:
        raise NotImplementedError

    def read(self, size: int | None = -1) -> bytes:
        chunk = self.source.read(size)
        if chunk:
            self.observe(chunk)
        return chunk


class MetricsRegistry:
    """Thread-safe registry to aggregate metrics."""

    def __init__(self):
        self._lock = threading.Lock()
        self.metrics: dict[str, dict[str, float]] = {}

    def update(self, name: str, size: int, duration: float):
        with self._lock:
            if name not in self.metrics:
                self.metrics[name] = {"bytes": 0, "time": 0.0}
            self.metrics[name]["bytes"] += size
            self.metrics[name]["time"] += duration

    def report(self):
        """Log the throughput per stage."""
        with self._lock:
            stats = []
            for name, data in self.metrics.items():
                mb = data["bytes"] / (1024 * 1024)
                if data["time"] > 0:
                    mb_s = mb / data["time"]
                    stats.append(f"{name}: {mb:.2f}MB in {data['time']:.2f}s ({mb_s:.2f} MB/s)")
            return " | ".join(stats)


class MeasuringStream(StreamWrapper):
    """Wraps a stream to measure read latency and throughput."""

    def __init__(self, source, name: str, registry: MetricsRegistry):
        super().__init__(source)
        self.name = name
        self.registry = registry

    def read(self, size: int | None = -1) -> bytes:
        start = time.perf_counter()
        chunk = self.source.read(size)
        duration = time.perf_counter() - start

        self.registry.update(self.name, len(chunk), duration)
        return chunk


class TqdmObserver(ObserverStream):
    def __init__(self, source: io.BufferedIOBase, pbar: Any | list[Any]):
        super().__init__(source)
        self.pbars = pbar if isinstance(pbar, list) else [pbar]
        self._lock = threading.Lock()

    def observe(self, chunk: bytes) -> None:
        n = len(chunk)
        with self._lock:
            for pbar in self.pbars:
                pbar.update(n)


class TransformStream(StreamWrapper, metaclass=ABCMeta):
    """
    Base class for streams that MODIFY data.
    """

    def __init__(self, source: io.BufferedIOBase):
        super().__init__(source)
        self._output_buffer = bytearray()
        self._eof = False

    @abstractmethod
    def _fill_buffer(self) -> bytes:
        """
        Read and transform data from the source stream.

        Subclasses must implement this method to read from the source, apply any transformation, and return the processed bytes.

        Example for a simple passthrough implementation::

            def _fill_buffer(self) -> bytes:
                return self.source.read(READ_CHUNK_SIZE)

        :returns: Transformed data chunk, or empty bytes when EOF is reached.
        """
        raise NotImplementedError

    def read(self, size: int | None = -1) -> bytes:
        target_size = size if size is not None else -1
        read_all = target_size == -1

        # Keep filling the output buffer until we have enough data to return or reach EOF
        while (read_all or len(self._output_buffer) < target_size) and (chunk := self._fill_buffer()):
            self._output_buffer.extend(chunk)

        # Return up to target_size bytes from the output buffer
        limit = min(len(self._output_buffer), target_size)
        ret = self._output_buffer[:limit]
        del self._output_buffer[:limit]
        return bytes(ret)


class Crypt4GHDecryptor(TransformStream):
    CIPHER_DIFF = crypt4gh.lib.CIPHER_DIFF
    CIPHER_SEGMENT_SIZE = crypt4gh.lib.CIPHER_SEGMENT_SIZE

    def __init__(self, source: io.BufferedIOBase, private_key: bytes):
        super().__init__(source)
        self._private_key = private_key
        self._session_keys: bytes | None = None
        self._header_parsed = False
        self._cipher_residue = bytearray()

    def _fill_buffer(self) -> bytes:
        if not self._header_parsed:
            self._read_header()

        return self._decrypt_next_segment()

    def _read_header(self) -> None:
        """Parse the Crypt4GH header to extract the session key for decryption."""
        try:
            keys = [(0, self._private_key, None)]
            # Decrypt header to extract session keys using crypt4gh library
            session_keys, edit_list = crypt4gh.header.deconstruct(self.source, keys, sender_pubkey=None)

            if edit_list is not None:
                raise ValueError("Edit lists in Crypt4GH headers are not supported!")

            if not session_keys:
                raise ValueError("No session keys found in Crypt4GH header")

            self._session_keys = session_keys
            self._header_parsed = True
        except Exception as e:
            raise OSError(f"Crypt4GH Header Error: {e}") from e

    def _decrypt_next_segment(self) -> bytes:
        """Decrypt and return the next segment, or None if EOF."""
        ciphersegment = self.source.read(self.CIPHER_SEGMENT_SIZE)
        if not ciphersegment:
            # EOF reached
            return b""
        elif len(ciphersegment) <= self.CIPHER_DIFF:
            # This means we have a truncated segment:
            # Even an empty plaintext segment would have CIPHER_DIFF bytes of overhead
            raise ValueError("Truncated cipher segment")
        else:
            return crypt4gh.lib.decrypt_block(ciphersegment, self._session_keys)


class Crypt4GHEncryptor(TransformStream):
    NONCE_LENGTH = 12
    SEGMENT_SIZE = crypt4gh.lib.SEGMENT_SIZE

    def __init__(self, source: io.BufferedIOBase, recipient_pubkey: bytes, sender_privkey: bytes | None = None):
        super().__init__(source)
        self._recipient_pubkey = recipient_pubkey
        self._sender_privkey = sender_privkey or bytes(PrivateKey.generate())
        self._session_key = os.urandom(32)
        self._header_sent = False
        self._plain_residue = bytearray()

    def _fill_buffer(self) -> bytes:
        if not self._header_sent:
            header = self._compose_header()
            self._header_sent = True

            return header

        return self._encrypt_next_segment()

    def _compose_header(self) -> bytes:
        keys = [(0, self._sender_privkey, self._recipient_pubkey)]

        header_content = crypt4gh.header.make_packet_data_enc(0, self._session_key)
        header_packets = crypt4gh.header.encrypt(header_content, keys)

        return crypt4gh.header.serialize(header_packets)

    def _encrypt_next_segment(self) -> bytes:
        segment = self.source.read(self.SEGMENT_SIZE)

        if len(segment) == 0:
            # EOF reached
            return b""

        nonce = os.urandom(self.NONCE_LENGTH)
        ciphertext = crypto_aead_chacha20poly1305_ietf_encrypt(segment, None, nonce, self._session_key)
        return nonce + ciphertext


class ValidatorObserver(ObserverStream):
    def __init__(self, source: io.BufferedIOBase, file_type: FileType, expected_checksum: str | None = None):
        super().__init__(source)
        self.file_type = file_type
        self.expected_checksum = expected_checksum
        self._hasher = hashlib.sha256()
        # if isa-l is available, use that:
        try:
            from isal import isal_zlib as izlib  # noqa: PLC0415

            self._decompressor = izlib.decompressobj(16 + zlib.MAX_WBITS)
        except ImportError:
            self._decompressor = zlib.decompressobj(16 + zlib.MAX_WBITS)
        self._fastq_line_buffer = b""
        self._fastq_line_count = 0
        self._bam_temp = None
        self._bam_path: str | None = None
        self._bytes_seen = 0

        if self.file_type == FileType.bam:
            self._bam_temp = tempfile.NamedTemporaryFile(suffix=".bam", delete=False)  # noqa: SIM115
            self._bam_path = self._bam_temp.name

    def observe(self, chunk: bytes) -> None:
        self._bytes_seen += len(chunk)
        self._hasher.update(chunk)

        if self.file_type == FileType.fastq:
            try:
                decompressed = self._decompressor.decompress(chunk)
                if decompressed:
                    self._validate_fastq_chunk(decompressed)
            except zlib.error as e:
                log.warning(f"Validator: GZIP error at byte {self._bytes_seen}: {e}")
        elif self.file_type == FileType.bam and self._bam_temp:
            self._bam_temp.write(chunk)

    def _validate_fastq_chunk(self, chunk: bytes) -> None:
        if self._fastq_line_buffer:
            chunk = self._fastq_line_buffer + chunk
            self._fastq_line_buffer = b""
        self._fastq_line_count += chunk.count(b"\n")
        r_idx = chunk.rfind(b"\n")
        if r_idx != -1:
            self._fastq_line_buffer = chunk[r_idx + 1 :]
        else:
            self._fastq_line_buffer = chunk

    @property
    def metrics(self) -> dict[str, Any]:
        stats: dict[str, Any] = {"checksum_sha256": self._hasher.hexdigest()}
        if self.file_type == FileType.fastq:
            stats["read_count"] = self._fastq_line_count // 4
        return stats

    def close(self) -> None:
        if not self.closed:
            if self.file_type == FileType.bam:
                if self._bam_temp:
                    with contextlib.suppress(Exception):
                        self._bam_temp.close()
                if self._bam_path and os.path.exists(self._bam_path):
                    with contextlib.suppress(OSError):
                        os.unlink(self._bam_path)
            super().close()

    def verify(self) -> None:
        if self.file_type == FileType.fastq:
            self._finalize_fastq()
        elif self.file_type == FileType.bam:
            self._finalize_bam()

        calculated = self._hasher.hexdigest()
        if self.expected_checksum and calculated != self.expected_checksum:
            raise ValueError(
                f"Checksum Mismatch! Seen {self._bytes_seen} bytes. Exp: {self.expected_checksum}, Got: {calculated}"
            )

    def _finalize_fastq(self) -> None:
        with contextlib.suppress(Exception):
            final = self._decompressor.flush()
            if final:
                self._validate_fastq_chunk(final)

        if self._fastq_line_buffer:
            self._fastq_line_count += 1

        if self._fastq_line_count > 0 and self._fastq_line_count % 4 != 0:
            raise ValueError(f"Invalid FASTQ: {self._fastq_line_count} lines (not divisible by 4)")

    def _finalize_bam(self) -> None:
        if self._bam_temp and self._bam_path:
            self._bam_temp.flush()
            try:
                pysam.AlignmentFile(self._bam_path, "rb", check_sq=False)
            except Exception as e:
                raise ValueError(f"BAM Validation Failed: {e}") from e


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
