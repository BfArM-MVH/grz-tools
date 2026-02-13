import contextlib
import hashlib
import io
import logging
import os
import tempfile
from typing import Any

import pysam
from grz_pydantic_models.submission.metadata import FileType

from . import ValidatingStream

try:
    from isal import isal_zlib as zlib
except ImportError:
    import zlib

log = logging.getLogger(__name__)


class ChecksumVerifier(ValidatingStream):
    def __init__(self, source: io.BufferedIOBase, algorithm: str = "sha256", expected_checksum: str | None = None):
        super().__init__(source)
        self.expected = expected_checksum
        self._hasher = hashlib.new(algorithm)
        self._bytes_seen = 0

    def observe(self, chunk: bytes) -> None:
        self._bytes_seen += len(chunk)
        self._hasher.update(chunk)

    def validate(self) -> None:
        calculated = self._hasher.hexdigest()
        if self.expected and calculated != self.expected:
            raise ValueError(f"Checksum mismatch! Exp: {self.expected}, Got: {calculated}")

    @property
    def metrics(self) -> dict[str, Any]:
        return {"checksum": self._hasher.hexdigest(), "size": self._bytes_seen}


class FastqVerifier(ValidatingStream):
    """
    Validates a decompressed FASTQ stream.
    """

    def __init__(self, source: io.BufferedIOBase):
        super().__init__(source)
        self._line_buffer = b""
        self._line_count = 0
        self._decompressor = zlib.decompressobj(16 + zlib.MAX_WBITS)

    def observe(self, chunk: bytes) -> None:
        try:
            decompressed = self._decompressor.decompress(chunk)
            if decompressed:
                self._count_lines(decompressed)
        except zlib.error as e:
            log.warning(f"FASTQ Validator: GZIP error: {e}")

    def _count_lines(self, chunk: bytes) -> None:
        if self._line_buffer:
            chunk = self._line_buffer + chunk
            self._line_buffer = b""

        self._line_count += chunk.count(b"\n")

        r_idx = chunk.rfind(b"\n")
        if r_idx != -1:
            self._line_buffer = chunk[r_idx + 1 :]
        else:
            self._line_buffer = chunk

    def validate(self) -> None:
        with contextlib.suppress(Exception):
            final = self._decompressor.flush()
            if final:
                self._count_lines(final)

        if self._line_buffer:
            self._line_count += 1

        if self._line_count > 0 and self._line_count % 4 != 0:
            raise ValueError(f"Invalid FASTQ: {self._line_count} lines (not divisible by 4)")

    @property
    def metrics(self) -> dict[str, Any]:
        return {"read_count": self._line_count // 4}


class BamVerifier(ValidatingStream):
    """
    Validates a BAM stream.
    Writes to a temp file because pysam requires random access.
    """

    def __init__(self, source: io.BufferedIOBase):
        super().__init__(source)
        self._temp = tempfile.NamedTemporaryFile(suffix=".bam", delete=False)
        self._path = self._temp.name

    def observe(self, chunk: bytes) -> None:
        self._temp.write(chunk)

    def validate(self) -> None:
        self._temp.close()
        try:
            pysam.AlignmentFile(self._path, "rb", check_sq=False)
        except Exception as e:
            raise ValueError(f"BAM Invalid: {e}") from e
        finally:
            if os.path.exists(self._path):
                with contextlib.suppress(OSError):
                    os.unlink(self._path)


class ValidatorObserver(ValidatingStream):
    """
    Composite observer that chains checksum and format validation.
    """

    def __init__(self, source: io.BufferedIOBase, file_type: FileType, expected_checksum: str | None = None):
        self._hasher = ChecksumVerifier(source, expected_checksum=expected_checksum)

        self._format_validator: ValidatingStream | None = None

        if file_type == FileType.fastq:
            self._format_validator = FastqVerifier(self._hasher)
        elif file_type == FileType.bam:
            self._format_validator = BamVerifier(self._hasher)

        upstream = self._format_validator if self._format_validator else self._hasher
        super().__init__(upstream)

    def observe(self, chunk: bytes) -> None:
        # No-op: The work is done by the upstream components (Hasher/Verifier)
        # when self.read() calls upstream.read()
        pass

    def validate(self) -> None:
        self._hasher.validate()

        if self._format_validator:
            self._format_validator.validate()

    @property
    def metrics(self) -> dict[str, Any]:
        data = self._hasher.metrics.copy()
        if self._format_validator:
            data.update(self._format_validator.metrics)
        return data
