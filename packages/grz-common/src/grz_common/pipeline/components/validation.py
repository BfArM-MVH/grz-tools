import contextlib
import hashlib
import io
import logging
import os
import tempfile
import zlib
from typing import Any

import pysam
from grz_pydantic_models.submission.metadata import FileType

from . import ObserverStream

log = logging.getLogger(__name__)


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
