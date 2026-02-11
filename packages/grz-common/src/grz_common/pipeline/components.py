"""Base classes and implementations for pipeline components."""

from __future__ import annotations

import hashlib
import io
import logging
import os
import tempfile
import threading
import zlib
from abc import ABCMeta, abstractmethod
from concurrent.futures import ThreadPoolExecutor
from typing import Any

import crypt4gh.header
import crypt4gh.lib
import pysam
from grz_pydantic_models.submission.metadata.v1 import FileType
from nacl.bindings import crypto_aead_chacha20poly1305_ietf_decrypt, crypto_aead_chacha20poly1305_ietf_encrypt
from nacl.public import PrivateKey

NONCE_LENGTH = 12
SEGMENT_SIZE = crypt4gh.lib.SEGMENT_SIZE
CIPHER_SEGMENT_SIZE = crypt4gh.lib.CIPHER_SEGMENT_SIZE
INT_WIDTH = 4
CRYPT4GH_MAGIC = b"crypt4gh"

log = logging.getLogger(__name__)


class PipelineError(Exception):
    """Base exception for pipeline errors."""

    def __init__(self, message: str, stage: str | None = None, cause: Exception | None = None):
        self.stage = stage
        self.cause = cause
        super().__init__(f"[{stage}] {message}" if stage else message)


class StreamWrapper(io.BufferedIOBase):
    """Base class wrapper."""

    def __init__(self, source: io.BufferedIOBase):
        self.source = source

    def readable(self) -> bool:
        return True

    def close(self):
        if not self.closed:
            self.source.close()
            super().close()


class TransformStream(StreamWrapper, metaclass=ABCMeta):
    """Base for streams that MODIFY data."""

    def __init__(self, source: io.BufferedIOBase):
        super().__init__(source)
        self._output_buffer = bytearray()
        self._eof = False

    @abstractmethod
    def _fill_buffer(self) -> None:
        raise NotImplementedError

    def read(self, size: int = -1) -> bytes:
        if size == -1:
            while not self._eof:
                self._fill_buffer()
            ret = bytes(self._output_buffer)
            self._output_buffer.clear()
            return ret

        if len(self._output_buffer) >= size:
            ret = self._output_buffer[:size]
            self._output_buffer = self._output_buffer[size:]
            return bytes(ret)

        while len(self._output_buffer) < size and not self._eof:
            self._fill_buffer()

        limit = min(len(self._output_buffer), size)
        ret = self._output_buffer[:limit]
        self._output_buffer = self._output_buffer[limit:]
        return bytes(ret)


class ObserverStream(StreamWrapper, metaclass=ABCMeta):
    """Base for streams that INSPECT data."""

    @abstractmethod
    def observe(self, chunk: bytes) -> None:
        raise NotImplementedError

    def read(self, size: int = -1) -> bytes:
        chunk = self.source.read(size)
        if chunk:
            self.observe(chunk)
        return chunk


class TqdmObserver(ObserverStream):
    def __init__(self, source: io.BufferedIOBase, pbar: Any):
        super().__init__(source)
        self.pbar = pbar
        self._lock = threading.Lock()

    def observe(self, chunk: bytes) -> None:
        with self._lock:
            self.pbar.update(len(chunk))


class Crypt4GHDecryptor(TransformStream):
    def __init__(self, source: io.BufferedIOBase, private_key: bytes):
        super().__init__(source)
        self._private_key = private_key
        self._session_key: bytes | None = None
        self._header_parsed = False
        self._header_buffer = bytearray()
        self._cipher_residue = bytearray()
        self._header_length = 0

    def _fill_buffer(self) -> None:
        # Phase 1: Parse Header
        if not self._header_parsed:
            chunk = self.source.read(io.DEFAULT_BUFFER_SIZE)
            if not chunk:
                if len(self._header_buffer) > 0:
                    log.error(f"Crypt4GH: Unexpected EOF. Header buffer has {len(self._header_buffer)} bytes.")
                    raise OSError("Unexpected EOF while reading Crypt4GH header")
                log.warning("Crypt4GH: Empty file detected (0 bytes read).")
                self._eof = True
                return

            self._header_buffer.extend(chunk)
            self._try_parse_header()
            return

        # Phase 2: Decrypt Body
        if len(self._cipher_residue) < CIPHER_SEGMENT_SIZE:
            chunk = self.source.read(io.DEFAULT_BUFFER_SIZE)
            if not chunk:
                self._eof = True
                self._flush_residue()
                return
            self._cipher_residue.extend(chunk)

        while len(self._cipher_residue) >= CIPHER_SEGMENT_SIZE:
            segment = self._cipher_residue[:CIPHER_SEGMENT_SIZE]
            self._cipher_residue = self._cipher_residue[CIPHER_SEGMENT_SIZE:]
            self._output_buffer.extend(self._decrypt_segment(segment))

    def _try_parse_header(self) -> None:
        if len(self._header_buffer) < 16:
            return

        try:
            with io.BytesIO(self._header_buffer) as f:
                magic = f.read(8)
                if magic != CRYPT4GH_MAGIC:
                    raise ValueError(f"Invalid Crypt4GH magic bytes: {magic!r}")

                _ = int.from_bytes(f.read(INT_WIDTH), "little")
                packet_count = int.from_bytes(f.read(INT_WIDTH), "little")

                for _ in range(packet_count):
                    if f.tell() + INT_WIDTH > len(self._header_buffer):
                        return
                    pkt_len = int.from_bytes(f.read(INT_WIDTH), "little")
                    if f.tell() + pkt_len - INT_WIDTH > len(self._header_buffer):
                        return
                    f.seek(pkt_len - INT_WIDTH, os.SEEK_CUR)

                self._header_length = f.tell()

            header_bytes = self._header_buffer[: self._header_length]
            keys = [(0, self._private_key, None)]

            session_keys, _ = crypt4gh.header.deconstruct(io.BytesIO(header_bytes), keys, sender_pubkey=None)

            self._session_key = session_keys[0]
            self._header_parsed = True

            log.debug(f"Crypt4GH Header Parsed. Residue size: {len(self._header_buffer) - self._header_length}")

            self._cipher_residue.extend(self._header_buffer[self._header_length :])
            self._header_buffer = None

        except Exception as e:
            raise OSError(f"Crypt4GH Header Error: {e}") from e

    def _decrypt_segment(self, segment: bytes) -> bytes:
        nonce = segment[:NONCE_LENGTH]
        return crypto_aead_chacha20poly1305_ietf_decrypt(segment[NONCE_LENGTH:], None, nonce, self._session_key)

    def _flush_residue(self) -> None:
        if self._cipher_residue:
            self._output_buffer.extend(self._decrypt_segment(bytes(self._cipher_residue)))
            self._cipher_residue.clear()


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

        chunk = self.source.read(io.DEFAULT_BUFFER_SIZE)
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
            self._output_buffer.extend(self._encrypt_segment(block))

    def _encrypt_segment(self, plaintext: bytes) -> bytes:
        nonce = os.urandom(NONCE_LENGTH)
        return nonce + crypto_aead_chacha20poly1305_ietf_encrypt(plaintext, None, nonce, self._session_key)


class ValidatorObserver(ObserverStream):
    def __init__(self, source: io.BufferedIOBase, file_type: FileType, expected_checksum: str | None = None):
        super().__init__(source)
        self.file_type = file_type
        self.expected_checksum = expected_checksum
        self._hasher = hashlib.sha256()
        self._decompressor = zlib.decompressobj(16 + zlib.MAX_WBITS)
        self._fastq_line_buffer = b""
        self._fastq_line_count = 0
        self._bam_temp = None
        self._bam_path: str | None = None

        if self.file_type == FileType.bam:
            self._bam_temp = tempfile.NamedTemporaryFile(suffix=".bam", delete=False)  # noqa: SIM115
            self._bam_path = self._bam_temp.name

    def observe(self, chunk: bytes) -> None:
        self._hasher.update(chunk)
        if self.file_type == FileType.fastq:
            try:
                decompressed = self._decompressor.decompress(chunk)
                if decompressed:
                    self._validate_fastq_chunk(decompressed)
            except zlib.error:
                # Suppress gzip errors mid-stream
                pass
        elif self.file_type == FileType.bam and self._bam_temp:
            self._bam_temp.write(chunk)

    def _validate_fastq_chunk(self, chunk: bytes) -> None:
        if self._fastq_line_buffer:
            chunk = self._fastq_line_buffer + chunk
            self._fastq_line_buffer = b""
        pos = 0
        while True:
            newline = chunk.find(b"\n", pos)
            if newline == -1:
                self._fastq_line_buffer = chunk[pos:]
                break
            self._fastq_line_count += 1
            pos = newline + 1

    @property
    def metrics(self) -> dict[str, Any]:
        stats: dict[str, Any] = {"checksum_sha256": self._hasher.hexdigest()}
        if self.file_type == FileType.fastq:
            stats["read_count"] = self._fastq_line_count // 4
        return stats

    def close(self) -> None:
        if not self.closed:
            if self.file_type == FileType.fastq:
                self._finalize_fastq()
            elif self.file_type == FileType.bam:
                self._finalize_bam()

            calculated = self._hasher.hexdigest()
            if self.expected_checksum and calculated != self.expected_checksum:
                raise ValueError(f"Checksum Mismatch! Exp: {self.expected_checksum}, Got: {calculated}")
            super().close()

    def _finalize_fastq(self) -> None:
        try:
            final = self._decompressor.flush()
            if final:
                self._validate_fastq_chunk(final)
        except Exception:
            pass
        if self._fastq_line_buffer:
            self._fastq_line_count += 1

        if self._fastq_line_count > 0 and self._fastq_line_count % 4 != 0:
            raise ValueError(f"Invalid FASTQ: {self._fastq_line_count} lines")

    def _finalize_bam(self) -> None:
        if self._bam_temp and self._bam_path:
            self._bam_temp.close()
            try:
                pysam.AlignmentFile(self._bam_path, "rb", check_sq=False)
            except Exception as e:
                raise ValueError(f"BAM Corrupt: {e}") from e
            finally:
                if os.path.exists(self._bam_path):
                    os.unlink(self._bam_path)


class S3Downloader(io.BufferedIOBase):
    def __init__(self, s3_client: Any, bucket: str, key: str):
        log.warning(f"S3Downloader Init: Bucket={bucket}, Key={key}")
        try:
            self.response = s3_client.get_object(Bucket=bucket, Key=key)
            self.stream = self.response["Body"]
            self.length = self.response.get("ContentLength", 0)
            log.warning(f"S3Downloader Opened: Length={self.length}")
        except Exception as e:
            log.error(f"S3Downloader Failed to Open {key}: {e}")
            raise

    def read(self, size: int = -1) -> bytes:
        try:
            data = self.stream.read(size) if size != -1 else self.stream.read()
            # if len(data) == 0:
            #     log.warning("S3Downloader: Read 0 bytes (EOF)")
            return data
        except Exception as e:
            log.error(f"S3Downloader Read Error: {e}")
            raise

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
        self._parts: list[dict] = []
        self._error: BaseException | None = None

    def upload(self, input_stream: io.BufferedIOBase) -> None:
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

    def _upload_part(self, uid: str, part_num: int, data: bytes) -> dict:
        local_md5 = hashlib.md5(data, usedforsecurity=False).digest()
        resp = self.s3.upload_part(Bucket=self.bucket, Key=self.key, UploadId=uid, PartNumber=part_num, Body=data)
        return {"PartNumber": part_num, "ETag": resp["ETag"], "local_md5": local_md5}

    def _collect_part(self, future) -> None:
        try:
            self._parts.append(future.result())
        except Exception as e:
            self._error = e

    def _calc_etag(self, parts: list) -> str:
        digests = [p["local_md5"] for p in parts if "local_md5" in p]
        if not digests:
            return ""
        combined = b"".join(digests)
        combined_hash = hashlib.md5(combined, usedforsecurity=False).hexdigest()
        return f"{combined_hash}-{len(digests)}"

    def abort(self) -> None:
        if self._upload_id:
            try:
                self.s3.abort_multipart_upload(Bucket=self.bucket, Key=self.key, UploadId=self._upload_id)
            except:
                pass
