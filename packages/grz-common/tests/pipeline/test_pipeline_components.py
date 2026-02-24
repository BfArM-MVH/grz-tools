"""Tests for the modular pipeline components."""

import gzip
import hashlib
from io import BytesIO

import pytest
from grz_common.pipeline.components import DataValidationError, ReadStream
from grz_common.pipeline.components.validation import ChecksumValidator, FastqValidator


class TestFastqValidator:
    """Tests for the FastqValidator stage."""

    def test_valid_fastq(self):
        """Test validation of a valid FASTQ file."""
        # Create valid FASTQ content (4 lines per record)
        fastq_content = b""
        for i in range(100):
            fastq_content += f"@read{i}\n".encode()
            fastq_content += b"ACGTACGTACGT\n"
            fastq_content += b"+\n"
            fastq_content += b"IIIIIIIIIIII\n"

        with (
            BytesIO(gzip.compress(fastq_content)) as f,
            ReadStream(f) as source,
            FastqValidator(mean_read_length_threshold=12) as validator,
        ):
            source >> validator
        metrics = validator.metrics

        assert metrics["line_count"] == 400  # 100 records * 4 lines
        assert metrics["read_count"] == 100
        assert metrics["mean_read_length"] == 12.0

    def test_invalid_fastq_line_count(self):
        """Test that non-multiple-of-4 line count is detected."""
        # Create invalid FASTQ (5 lines instead of 4)
        fastq_content = b"@read1\nACGT\n+\nIIII\nextra_line\n"

        with pytest.raises(DataValidationError, match=r"Invalid FASTQ"):
            with (
                BytesIO(gzip.compress(fastq_content)) as f,
                ReadStream(f) as source,
                FastqValidator(mean_read_length_threshold=12) as validator,
            ):
                source >> validator

    def test_invalid_fastq_seq_qual_length_mismatch(self):
        """Test that sequence and quality length match."""
        fastq_content = b"@read1\nACGT\n+\nII\n"

        with pytest.raises(DataValidationError, match=r"Invalid FASTQ"):
            with (
                BytesIO(gzip.compress(fastq_content)) as f,
                ReadStream(f) as source,
                FastqValidator(mean_read_length_threshold=12) as validator,
            ):
                source >> validator


class TestRawChecksumValidator:
    """Tests for the ChecksumValidator."""

    def test_valid_checksum(self):
        """Test validation with correct checksum."""
        data = b"Test data for checksum validation"
        expected_checksum = hashlib.sha256(data).hexdigest()

        with (
            BytesIO(data) as f,
            ReadStream(f) as source,
            ChecksumValidator(expected_checksum=expected_checksum) as validator,
        ):
            source >> validator

    def test_invalid_checksum(self):
        """Test that checksum mismatch is detected."""
        data = b"Test data for checksum validation"

        with pytest.raises(DataValidationError, match=r"Checksum mismatch!"):
            with (
                BytesIO(data) as f,
                ReadStream(f) as source,
                ChecksumValidator(expected_checksum="0" * 64) as validator,
            ):
                source >> validator
