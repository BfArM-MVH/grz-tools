"""Tests for the modular pipeline components."""

import gzip
import hashlib
import io

from grz_common.pipeline.base import PipelineContext
from grz_common.pipeline.compressors import GzipDecompressor
from grz_common.pipeline.validators import FastqValidator, RawChecksumValidator


class TestGzipDecompressor:
    """Tests for the GzipDecompressor stage."""

    def test_decompress_gzipped_data(self):
        """Test decompression of gzipped data."""
        original = b"Hello, World! " * 1000
        compressed = gzip.compress(original)

        decompressor = GzipDecompressor()
        context = PipelineContext()
        decompressor.initialize(context)

        # Process in chunks
        chunk_size = 100
        result = io.BytesIO()
        for i in range(0, len(compressed), chunk_size):
            chunk = compressed[i : i + chunk_size]
            decompressed = decompressor.process(chunk)
            result.write(decompressed)

        # Flush remaining
        result.write(decompressor.flush())
        decompressor.finalize()

        assert result.getvalue() == original
        assert decompressor.is_gzipped is True

    def test_passthrough_non_gzipped_data(self):
        """Test that non-gzipped data passes through unchanged."""
        original = b"This is not gzipped data"

        decompressor = GzipDecompressor()
        context = PipelineContext()
        decompressor.initialize(context)

        result = decompressor.process(original)
        result += decompressor.flush()
        decompressor.finalize()

        assert result == original
        assert decompressor.is_gzipped is False


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

        validator = FastqValidator(mean_read_length_threshold=0)
        context = PipelineContext()
        validator.initialize(context)

        # Feed in chunks
        chunk_size = 50
        for i in range(0, len(fastq_content), chunk_size):
            validator.observe(fastq_content[i : i + chunk_size])

        validator.finalize()

        assert not context.has_errors()
        assert validator.line_count == 400  # 100 records * 4 lines
        assert validator.read_count == 100
        assert validator.mean_read_length == 12.0

    def test_invalid_fastq_line_count(self):
        """Test that non-multiple-of-4 line count is detected."""
        # Create invalid FASTQ (5 lines instead of 4)
        fastq_content = b"@read1\nACGT\n+\nIIII\nextra_line\n"

        validator = FastqValidator(mean_read_length_threshold=0)
        context = PipelineContext()
        validator.initialize(context)
        validator.observe(fastq_content)
        validator.finalize()

        assert context.has_errors()
        assert any("not a multiple of 4" in e for e in context.errors)

    def test_fastq_with_gzip_decompressor(self):
        """Test FASTQ validation with upstream decompressor."""
        # Create valid FASTQ content
        fastq_content = b""
        for i in range(50):
            fastq_content += f"@read{i}\n".encode()
            fastq_content += b"ACGTACGTACGTACGT\n"
            fastq_content += b"+\n"
            fastq_content += b"IIIIIIIIIIIIIIII\n"

        # Gzip compress
        compressed = gzip.compress(fastq_content)

        # Setup pipeline stages
        context = PipelineContext()
        decompressor = GzipDecompressor()
        validator = FastqValidator(mean_read_length_threshold=0)

        decompressor.initialize(context)
        validator.initialize(context)

        # Process through decompressor -> validator
        chunk_size = 100
        for i in range(0, len(compressed), chunk_size):
            chunk = compressed[i : i + chunk_size]
            decompressed = decompressor.process(chunk)
            if decompressed:
                validator.observe(decompressed)

        # Flush
        final_decompressed = decompressor.flush()
        if final_decompressed:
            validator.observe(final_decompressed)

        decompressor.finalize()
        validator.finalize()

        assert not context.has_errors()
        assert validator.line_count == 200  # 50 records * 4 lines
        assert validator.read_count == 50
        assert validator.mean_read_length == 16.0


class TestRawChecksumValidator:
    """Tests for the RawChecksumValidator stage."""

    def test_valid_checksum(self):
        """Test validation with correct checksum."""
        data = b"Test data for checksum validation"
        expected_checksum = hashlib.sha256(data).hexdigest()
        expected_size = len(data)

        validator = RawChecksumValidator(
            expected_checksum=expected_checksum,
            expected_size=expected_size,
        )
        context = PipelineContext()
        validator.initialize(context)
        validator.observe(data)
        validator.finalize()

        assert not context.has_errors()
        assert validator.calculated_checksum == expected_checksum

    def test_invalid_checksum(self):
        """Test that checksum mismatch is detected."""
        data = b"Test data for checksum validation"

        validator = RawChecksumValidator(
            expected_checksum="0" * 64,  # Wrong checksum
            expected_size=len(data),
        )
        context = PipelineContext()
        validator.initialize(context)
        validator.observe(data)
        validator.finalize()

        assert context.has_errors()
        assert any("Checksum mismatch" in e for e in context.errors)

    def test_invalid_size(self):
        """Test that size mismatch is detected."""
        data = b"Test data"

        validator = RawChecksumValidator(
            expected_checksum=hashlib.sha256(data).hexdigest(),
            expected_size=len(data) + 100,  # Wrong size
        )
        context = PipelineContext()
        validator.initialize(context)
        validator.observe(data)
        validator.finalize()

        assert context.has_errors()
        assert any("size mismatch" in e for e in context.errors)
