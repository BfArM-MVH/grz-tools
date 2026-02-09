"""Validation pipeline stages."""

from __future__ import annotations

import contextlib
import hashlib
import tempfile
from pathlib import Path
from typing import Any

from .base import StreamObserver


class RawChecksumValidator(StreamObserver):
    """Validates checksum and file size of raw data."""

    def __init__(
        self,
        expected_checksum: str | None = None,
        expected_size: int | None = None,
        checksum_type: str = "sha256",
        name: str | None = None,
    ):
        """
        Initialize the checksum validator.

        :param expected_checksum: Expected checksum (hex string)
        :param expected_size: Expected file size in bytes
        :param checksum_type: Type of checksum (only sha256 supported)
        :param name: Stage name for logging
        """
        super().__init__(name or "RawChecksumValidator")
        self._expected_checksum = expected_checksum
        self._expected_size = expected_size
        self._checksum_type = checksum_type.lower()

        if self._checksum_type != "sha256":
            raise ValueError(f"Unsupported checksum type: {checksum_type}")

        self._hasher = hashlib.sha256()
        self._bytes_observed = 0

    def observe(self, data: bytes) -> None:
        """Update checksum with observed data."""
        self._hasher.update(data)
        self._bytes_observed += len(data)

    def finalize(self) -> None:
        """Validate checksum and size, adding errors to context if invalid."""
        calculated = self._hasher.hexdigest()

        # store calculated values in context
        self.context.checksums["sha256"] = calculated
        self.context.metadata["validated_size"] = self._bytes_observed

        # validate checksum
        if self._expected_checksum and calculated != self._expected_checksum:
            error = f"Checksum mismatch: expected {self._expected_checksum}, got {calculated}"
            self.context.add_error(error)
            self._log.error(error)

        # validate size
        if self._expected_size is not None and self._bytes_observed != self._expected_size:
            error = f"File size mismatch: expected {self._expected_size} bytes, got {self._bytes_observed} bytes"
            self.context.add_error(error)
            self._log.error(error)

    @property
    def calculated_checksum(self) -> str:
        """Return the calculated checksum."""
        return self._hasher.hexdigest()

    @property
    def bytes_observed(self) -> int:
        """Return bytes observed."""
        return self._bytes_observed

    def get_result(self) -> dict[str, Any]:
        """Return validation results."""
        return {
            "checksum": self._hasher.hexdigest(),
            "size": self._bytes_observed,
            "valid": not self.context.has_errors(),
        }


class FastqValidator(StreamObserver):
    """Streaming FASTQ format validator for decompressed data."""

    def __init__(
        self,
        mean_read_length_threshold: int = 0,
        name: str | None = None,
    ):
        """
        Initialize the FASTQ validator.

        :param mean_read_length_threshold: Minimum mean read length (exclusive)
        :param name: Stage name for logging
        """
        super().__init__(name or "FastqValidator")
        self._mean_read_length_threshold = mean_read_length_threshold

        # validation state
        self._line_buffer = b""
        self._line_number = 0
        self._total_read_length = 0
        self._total_reads = 0
        self._bytes_observed = 0

    def observe(self, data: bytes) -> None:
        """
        Process observed DECOMPRESSED data for FASTQ validation.

        :param data: Decompressed FASTQ data bytes
        """
        self._bytes_observed += len(data)

        # process lines directly (data should already be decompressed)
        self._validate_chunk(data)

    def _validate_chunk(self, data: bytes) -> None:
        """Validate a chunk of decompressed FASTQ data."""
        # combine with leftover from previous chunk
        if self._line_buffer:
            data = self._line_buffer + data
            self._line_buffer = b""

        # find line boundaries without splitting into many objects
        # this is significantly faster than data.split(b"\n")
        pos = 0
        data_len = len(data)
        line_number = self._line_number
        total_read_length = self._total_read_length
        total_reads = self._total_reads

        while pos < data_len:
            # find next newline
            newline_pos = data.find(b"\n", pos)

            if newline_pos == -1:
                # no more newlines, save remainder for next chunk
                self._line_buffer = data[pos:]
                break

            # sequence lines are every 4th line starting from index 1 (0-indexed)
            if (line_number % 4) == 1:
                # calculate line length without carriage return
                line_end = newline_pos
                if line_end > pos and data[line_end - 1 : line_end] == b"\r":
                    line_end -= 1
                total_read_length += line_end - pos
                total_reads += 1

            line_number += 1
            pos = newline_pos + 1

        self._line_number = line_number
        self._total_read_length = total_read_length
        self._total_reads = total_reads

    def finalize(self) -> None:
        """Finalize validation and add any errors to context."""
        # process remaining buffer
        if self._line_buffer:
            line_text = self._line_buffer.decode("utf-8", errors="replace").rstrip("\r")
            if line_text:
                if (self._line_number % 4) == 1:
                    self._total_read_length += len(line_text)
                    self._total_reads += 1
                self._line_number += 1

        # store stats in context
        self.context.metadata["fastq_line_count"] = self._line_number
        self.context.metadata["fastq_read_count"] = self._total_reads

        mean_length = self.mean_read_length  # use property to avoid duplication
        if self._total_reads > 0:
            self.context.metadata["fastq_mean_read_length"] = mean_length

        # validate line count
        if self._line_number > 0 and self._line_number % 4 != 0:
            error = f"Number of lines is not a multiple of 4! Found {self._line_number} lines."
            self.context.add_error(error)
            self._log.error(error)

        # validate mean read length
        if self._total_reads > 0 and mean_length <= self._mean_read_length_threshold:
            error = f"Mean read length must be > {self._mean_read_length_threshold}bp, calculated {mean_length:.2f}."
            self.context.add_error(error)
            self._log.error(error)

        self._log.debug(f"FASTQ validation complete: {self._line_number} lines, {self._total_reads} reads")

    @property
    def line_count(self) -> int:
        """Return number of lines processed."""
        return self._line_number

    @property
    def read_count(self) -> int:
        """Return number of reads processed."""
        return self._total_reads

    @property
    def mean_read_length(self) -> float:
        """Return mean read length."""
        if self._total_reads == 0:
            return 0.0
        return self._total_read_length / self._total_reads

    def get_result(self) -> dict[str, Any]:
        """Return validation results."""
        return {
            "line_count": self._line_number,
            "read_count": self._total_reads,
            "mean_read_length": self.mean_read_length,
            "valid": not self.context.has_errors(),
        }


class BamValidator(StreamObserver):
    """
    Streaming BAM validator.

    Due to pysam's requirement for file paths, this validator writes data
    to a temporary file and validates when finalize() is called.

    Checks:
    - Header presence (warning only, not an error)
    - Secondary alignments (warning only)
    - Hard-clipped bases (warning only)

    Note: Current BAM validation only produces warnings, not errors.

    Usage:
        validator = BamValidator()
        validator.initialize(context)
        validator.observe(bam_data)
        validator.finalize()  # Performs validation
    """

    def __init__(self, name: str | None = None):
        """Initialize the BAM validator."""
        super().__init__(name or "BamValidator")
        self._temp_file: Any = None
        self._temp_path: str | None = None
        self._bytes_written = 0

    def initialize(self, context: Any) -> None:
        """Create temporary file for BAM data."""
        super().initialize(context)
        self._temp_file = tempfile.NamedTemporaryFile(suffix=".bam", delete=False)  # noqa: SIM115
        self._temp_path = self._temp_file.name

    def observe(self, data: bytes) -> None:
        """Write data to temporary file."""
        if self._temp_file is not None:
            self._temp_file.write(data)
            self._bytes_written += len(data)

    def finalize(self) -> None:
        """Validate the BAM file."""
        if self._temp_file is None or self._temp_path is None:
            return

        # close temp file for reading
        self._temp_file.close()

        try:
            import pysam  # noqa: PLC0415

            # disable SQ header enforcement (not needed for uBAM)
            bam_file = pysam.AlignmentFile(self._temp_path, mode="rb", check_sq=False)

            header = bam_file.header.to_dict()
            # HD key is usually present and safe
            concerning_keys = header.keys() - {"HD"}
            if concerning_keys:
                self._log.warning("Detected a header in BAM, ensure it contains no private information!")

            secondary_warned = False
            hard_clipped_warned = False

            # need until_eof since BAM won't have an index
            for read in bam_file.fetch(until_eof=True):
                if not secondary_warned and read.is_secondary:
                    self._log.warning(
                        "Detected secondary alignment in BAM. Consider filtering to save bandwidth and storage."
                    )
                    secondary_warned = True

                if not hard_clipped_warned and not read.is_secondary:
                    hard_clipped_count = read.get_cigar_stats()[0][5]
                    if hard_clipped_count:
                        self._log.warning(
                            "Detected hard-clipped bases in primary alignment. "
                            "This is a loss of information from raw reads."
                        )
                        hard_clipped_warned = True

            bam_file.close()

            self.context.metadata["bam_bytes_validated"] = self._bytes_written
            self._log.debug(f"BAM validation complete ({self._bytes_written} bytes)")

        except ImportError:
            self._log.warning("pysam not available, skipping BAM validation")
        except Exception as e:
            error = f"BAM validation failed: {e}"
            self.context.add_error(error)
            self._log.error(error)
        finally:
            # clean up temp file
            self._cleanup_temp()

    def abort(self) -> None:
        """Clean up on abort."""
        if self._temp_file is not None:
            with contextlib.suppress(Exception):
                self._temp_file.close()
        self._cleanup_temp()

    def _cleanup_temp(self) -> None:
        """Remove temporary file."""
        if self._temp_path:
            with contextlib.suppress(Exception):
                Path(self._temp_path).unlink(missing_ok=True)
            self._temp_path = None

    def get_result(self) -> dict[str, Any]:
        """Return validation results."""
        return {
            "bytes_validated": self._bytes_written,
            "valid": not self.context.has_errors(),
        }
