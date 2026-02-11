"""Base classes and implementations for pipeline components."""

from __future__ import annotations

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
from nacl.bindings import (
    crypto_aead_chacha20poly1305_ietf_decrypt,
    crypto_aead_chacha20poly1305_ietf_encrypt,
)
from nacl.public import PrivateKey

log = logging.getLogger(__name__)

NONCE_LENGTH = 12
SEGMENT_SIZE = crypt4gh.lib.SEGMENT_SIZE
CIPHER_SEGMENT_SIZE = crypt4gh.lib.CIPHER_SEGMENT_SIZE
INT_WIDTH = 4
CRYPT4GH_MAGIC = b"crypt4gh"

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


class TransformStream(StreamWrapper, metaclass=ABCMeta):
    """
    Base class for streams that MODIFY data.
    """

    def __init__(self, source: io.BufferedIOBase):
        super().__init__(source)
        self._output_buffer = bytearray()
        self._eof = False

    @abstractmethod
    def _fill_buffer(self) -> None:
        raise NotImplementedError

    def read(self, size: int | None = -1) -> bytes:
        target_size = size if size is not None else -1

        if target_size == -1:
            while not self._eof:
                self._fill_buffer()
            ret = bytes(self._output_buffer)
            self._output_buffer.clear()
            return ret

        if len(self._output_buffer) >= target_size:
            ret = self._output_buffer[:target_size]
            self._output_buffer = self._output_buffer[target_size:]
            return bytes(ret)

        while len(self._output_buffer) < target_size and not self._eof:
            self._fill_buffer()

        limit = min(len(self._output_buffer), target_size)
        ret = self._output_buffer[:limit]
        self._output_buffer = self._output_buffer[limit:]
        return bytes(ret)


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


class Crypt4GHDecryptor(TransformStream):
    def __init__(self, source: io.BufferedIOBase, private_key: bytes):
        super().__init__(source)
        self._private_key = private_key
        self._session_key: bytes | None = None
        self._header_parsed = False
        self._header_buffer = bytearray()
        self._cipher_residue = bytearray()

    def _fill_buffer(self) -> None:
        if not self._header_parsed:
            self._process_header()
            if not self._header_parsed:
                return

        self._process_body()

    def _process_header(self) -> None:
        required_min = 16
        while len(self._header_buffer) < required_min:
            chunk = self.source.read(required_min - len(self._header_buffer))
            if not chunk:
                if len(self._header_buffer) > 0:
                    raise OSError(f"Crypt4GH: Unexpected EOF reading header (got {len(self._header_buffer)} bytes)")
                self._eof = True
                return
            self._header_buffer.extend(chunk)

        self._try_parse_header()

    def _process_body(self) -> None:
        while len(self._cipher_residue) < CIPHER_SEGMENT_SIZE:
            chunk = self.source.read(READ_CHUNK_SIZE)
            if not chunk:
                self._eof = True
                break
            self._cipher_residue.extend(chunk)

        while len(self._cipher_residue) >= CIPHER_SEGMENT_SIZE:
            segment = self._cipher_residue[:CIPHER_SEGMENT_SIZE]
            self._cipher_residue = self._cipher_residue[CIPHER_SEGMENT_SIZE:]
            self._output_buffer.extend(self._decrypt_segment(segment))

        if self._eof and self._cipher_residue:
            self._output_buffer.extend(self._decrypt_segment(bytes(self._cipher_residue)))
            self._cipher_residue.clear()

    def _try_parse_header(self) -> None:
        try:
            with io.BytesIO(self._header_buffer) as f:
                if f.read(8) != CRYPT4GH_MAGIC:
                    raise ValueError("Invalid Crypt4GH Magic")

                _ = int.from_bytes(f.read(INT_WIDTH), "little")
                packet_count = int.from_bytes(f.read(INT_WIDTH), "little")

                for _ in range(packet_count):
                    if f.tell() + INT_WIDTH > len(self._header_buffer):
                        self._read_more_header(INT_WIDTH)
                        return

                    f_pos = f.tell()
                    pkt_len = int.from_bytes(self._header_buffer[f_pos : f_pos + INT_WIDTH], "little")
                    f.seek(INT_WIDTH, os.SEEK_CUR)

                    needed = pkt_len - INT_WIDTH
                    if f.tell() + needed > len(self._header_buffer):
                        self._read_more_header(needed)
                        return

                    f.seek(needed, os.SEEK_CUR)

                header_length = f.tell()

            header_bytes = self._header_buffer[:header_length]
            keys = [(0, self._private_key, None)]
            session_keys, _ = crypt4gh.header.deconstruct(io.BytesIO(header_bytes), keys, sender_pubkey=None)

            if not session_keys:
                raise ValueError("No session keys found in Crypt4GH header")

            self._session_key = session_keys[0]
            self._header_parsed = True

            self._cipher_residue.extend(self._header_buffer[header_length:])
            self._header_buffer.clear()

        except Exception as e:
            raise OSError(f"Crypt4GH Header Error: {e}") from e

    def _read_more_header(self, size: int) -> None:
        chunk = self.source.read(size)
        if chunk:
            self._header_buffer.extend(chunk)

    def _decrypt_segment(self, segment: bytes) -> bytes:
        assert self._session_key is not None, "Session key not initialized"
        nonce = segment[:NONCE_LENGTH]
        ciphertext = segment[NONCE_LENGTH:]
        return crypto_aead_chacha20poly1305_ietf_decrypt(bytes(ciphertext), None, bytes(nonce), self._session_key)


class Crypt4GHEncryptor(TransformStream):
    def __init__(self, source: io.BufferedIOBase, recipient_pubkey: bytes, sender_privkey: bytes | None = None):
        super().__init__(source)
        self._recipient_pubkey = recipient_pubkey
        self._sender_privkey = sender_privkey or bytes(PrivateKey.generate())
        self._session_key = os.urandom(32)
        self._header_sent = False
        self._plain_residue = bytearray()

    def _fill_buffer(self) -> None:
        if not self._header_sent:
            keys = [(0, self._sender_privkey, self._recipient_pubkey)]
            header_content = crypt4gh.header.make_packet_data_enc(0, self._session_key)
            header_packets = crypt4gh.header.encrypt(header_content, keys)
            self._output_buffer.extend(crypt4gh.header.serialize(header_packets))
            self._header_sent = True

        chunk = self.source.read(READ_CHUNK_SIZE)
        if not chunk:
            self._eof = True
            if self._plain_residue:
                self._output_buffer.extend(self._encrypt_segment(bytes(self._plain_residue)))
                self._plain_residue.clear()
            return

        self._plain_residue.extend(chunk)
        while len(self._plain_residue) >= SEGMENT_SIZE:
            block = self._plain_residue[:SEGMENT_SIZE]
            self._plain_residue = self._plain_residue[SEGMENT_SIZE:]
            self._output_buffer.extend(self._encrypt_segment(bytes(block)))

    def _encrypt_segment(self, plaintext: bytes) -> bytes:
        nonce = os.urandom(NONCE_LENGTH)
        ciphertext = crypto_aead_chacha20poly1305_ietf_encrypt(bytes(plaintext), None, nonce, self._session_key)
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


class S3MultipartUploader:
    def __init__(  # noqa: PLR0913
        self,
        s3_client: Any,
        bucket: str,
        key: str,
        part_size: int | None = None,
        max_threads: int = 4,
    ):
        self.s3 = s3_client
        self.bucket = bucket
        self.key = key
        self.part_size = part_size
        self.executor = ThreadPoolExecutor(max_workers=max_threads)
        self._upload_id: str | None = None
        self._parts: list[dict[str, Any]] = []
        self._error: BaseException | None = None

    def upload(self, input_stream: io.BufferedIOBase) -> None:
        log.info(f"S3Uploader: Starting upload to s3://{self.bucket}/{self.key} (Part size: {self.part_size})")
        try:
            resp = self.s3.create_multipart_upload(Bucket=self.bucket, Key=self.key)
            self._upload_id = resp["UploadId"]

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

                while len(futures) > 4:
                    self._collect_part(futures.pop(0))

            for f in futures:
                self._collect_part(f)

            if self._error:
                raise self._error

            self._parts.sort(key=lambda x: x["PartNumber"])
            expected = self._calc_etag(self._parts)

            complete = self.s3.complete_multipart_upload(
                Bucket=self.bucket,
                Key=self.key,
                UploadId=self._upload_id,
                MultipartUpload={"Parts": [{"PartNumber": p["PartNumber"], "ETag": p["ETag"]} for p in self._parts]},
            )

            server_etag = complete.get("ETag", "").strip('"')
            if expected and server_etag != expected:
                raise OSError(f"ETag mismatch! Exp: {expected}, Got: {server_etag}")

        except Exception as e:
            self.abort()
            raise e
        finally:
            self.executor.shutdown(wait=False)

    def _upload_part(self, uid: str, part_num: int, data: bytes) -> dict[str, Any]:
        local_md5 = hashlib.md5(data, usedforsecurity=False).digest()
        resp = self.s3.upload_part(Bucket=self.bucket, Key=self.key, UploadId=uid, PartNumber=part_num, Body=data)
        return {"PartNumber": part_num, "ETag": resp["ETag"], "local_md5": local_md5}

    def _collect_part(self, future) -> None:
        try:
            self._parts.append(future.result())
        except Exception as e:
            self._error = e

    def _calc_etag(self, parts: list[dict[str, Any]]) -> str:
        digests = [p["local_md5"] for p in parts if "local_md5" in p]
        if not digests:
            return ""
        combined = b"".join(digests)
        combined_hash = hashlib.md5(combined, usedforsecurity=False).hexdigest()
        return f"{combined_hash}-{len(digests)}"

    def abort(self) -> None:
        if self._upload_id:
            with contextlib.suppress(Exception):
                self.s3.abort_multipart_upload(Bucket=self.bucket, Key=self.key, UploadId=self._upload_id)
