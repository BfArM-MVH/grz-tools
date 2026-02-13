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

try:
    import numpy as np
    from numba import jit

    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False
    np = None
    jit = None

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


class FastqValidator(ValidatingStream):
    """
    Validates a decompressed FASTQ stream.
    """

    def __init__(self, source: io.BufferedIOBase, mean_read_length_threshold: float | None = None):
        super().__init__(source)
        self._decompressor = zlib.decompressobj(16 + zlib.MAX_WBITS)
        self._processor = self._process_chunk_numba if HAS_NUMBA else self._process_chunk_python
        self._threshold = mean_read_length_threshold

        self._line_state = 0
        self._total_read_len = 0
        self._read_count = 0
        self._current_seq_len = 0
        self._current_line_len = 0
        self._last_chunk_byte = 10  # b'\n'

        if HAS_NUMBA:
            log.debug("FastqValidator using Numba optimization.")

    def observe(self, chunk: bytes) -> None:
        try:
            decompressed = self._decompressor.decompress(chunk)
            if decompressed:
                self._processor(decompressed)
        except zlib.error as e:
            log.warning(f"FASTQ Validator: GZIP error: {e}")

    def _process_chunk_numba(self, chunk: bytes) -> None:
        chunk_arr = np.frombuffer(chunk, dtype=np.uint8)

        (
            status,
            self._line_state,
            self._read_count,
            self._total_read_len,
            self._current_seq_len,
            self._current_line_len,
            self._last_chunk_byte,
        ) = _validate_chunk_numba(
            chunk_arr,
            self._line_state,
            self._read_count,
            self._total_read_len,
            self._current_seq_len,
            self._current_line_len,
            self._last_chunk_byte,
        )

        if status == 1:
            raise ValueError(
                f"FASTQ Record {self._read_count + 1} Invalid: "
                f"Seq len ({self._current_seq_len}) != Qual len (inconsistent)"
            )

    def _process_chunk_python(self, chunk: bytes) -> None:  # noqa: C901
        start = 0
        n_len = len(chunk)

        while True:
            pos = chunk.find(b"\n", start)

            if pos == -1:
                self._current_line_len += n_len - start
                if n_len > 0:
                    self._last_chunk_byte = chunk[-1]
                break

            self._current_line_len += pos - start

            # Check for \r
            is_cr = False
            if pos > 0:
                if chunk[pos - 1] == 13:  # b'\r'
                    is_cr = True
            elif self._last_chunk_byte == 13:
                is_cr = True

            final_len = self._current_line_len
            if is_cr:
                final_len -= 1

            final_len = max(0, final_len)

            if self._line_state == 1:
                self._current_seq_len = final_len
                self._total_read_len += final_len
                self._line_state = 2

            elif self._line_state == 3:
                if final_len != self._current_seq_len:
                    raise ValueError(
                        f"FASTQ Record {self._read_count + 1} Invalid: "
                        f"Seq len ({self._current_seq_len}) != Qual len ({final_len})"
                    )
                self._read_count += 1
                self._line_state = 0

            else:
                self._line_state += 1

            self._current_line_len = 0
            self._last_chunk_byte = 10
            start = pos + 1

    def validate(self) -> None:
        with contextlib.suppress(Exception):
            final = self._decompressor.flush()
            if final:
                if HAS_NUMBA:
                    self._process_chunk_numba(final)
                else:
                    self._process_chunk_python(final)

        if self._current_line_len > 0:
            raise ValueError("Invalid FASTQ: Unexpected EOF (incomplete line). File may be truncated.")

        if self._line_state != 0:
            raise ValueError(
                f"Invalid FASTQ: Unexpected EOF (incomplete record). Finished at state {self._line_state}."
            )

        if self._threshold is not None and self._read_count > 0:
            mean_length = self._total_read_len / self._read_count
            if mean_length < self._threshold:
                raise ValueError(f"Mean read length ({mean_length:.2f}) is below threshold ({self._threshold})")

    @property
    def metrics(self) -> dict[str, Any]:
        return {"read_count": self._read_count, "total_bases": self._total_read_len}


if HAS_NUMBA:

    @jit(nopython=True, nogil=True, cache=True)
    def _validate_chunk_numba(chunk_data, line_state, read_count, total_len, curr_seq_len, curr_line_len, last_byte):  # noqa: C901, PLR0913
        """
        Validates a single chunk of FASTQ data using Numba.

        :return: (status [0 on success, 1 on error], line_state, read_count, total_len, curr_seq_len, curr_line_len, last_byte)
        """
        n = len(chunk_data)

        for i in range(n):
            byte = chunk_data[i]

            if byte != 10:  # not b'\n'
                curr_line_len += 1
                continue

            # otherwise found newline
            final_len = curr_line_len

            # handle \r\n style newlines
            if i > 0:
                if chunk_data[i - 1] == 13:  # \r
                    final_len -= 1
            elif last_byte == 13:
                final_len -= 1

            # Clamp length to 0
            final_len = max(final_len, 0)

            if line_state == 0:
                line_state = 1
            elif line_state == 1:
                curr_seq_len = final_len
                total_len += final_len
                line_state = 2
            elif line_state == 2:
                line_state = 3
            elif line_state == 3:
                if final_len != curr_seq_len:
                    # quality length does not match sequence length
                    return 1, line_state, read_count, total_len, curr_seq_len, curr_line_len, last_byte

                read_count += 1
                line_state = 0

            curr_line_len = 0

        # update last_byte for next chunk
        if n > 0:
            last_byte = chunk_data[n - 1]

        return 0, line_state, read_count, total_len, curr_seq_len, curr_line_len, last_byte


class BamValidator(ValidatingStream):
    """
    Validates a BAM stream.
    Writes to a temp file because pysam requires random access.
    """

    def __init__(self, source: io.BufferedIOBase):
        super().__init__(source)
        self._temp = tempfile.NamedTemporaryFile(suffix=".bam", delete=False)  # noqa: SIM115
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

    def __init__(self, source: io.BufferedIOBase, file_type: FileType, expected_checksum: str | None = None, **kwargs):
        self._hasher = ChecksumVerifier(source, expected_checksum=expected_checksum)
        self._format_validator: ValidatingStream | None = None

        if file_type == FileType.fastq:
            threshold = kwargs.get("mean_read_length_threshold")
            self._format_validator = FastqValidator(self._hasher, mean_read_length_threshold=threshold)
        elif file_type == FileType.bam:
            self._format_validator = BamValidator(self._hasher)

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
