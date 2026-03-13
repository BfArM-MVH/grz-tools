"""Tests for the modular pipeline components."""

import gzip
import hashlib
from io import BytesIO

import pytest
from grz_common.pipeline.components import DataValidationError, Observer, ReadStream, Transformer
from grz_common.pipeline.components.validation import ChecksumValidator, FastqValidator


class _CountingChunkTransformer(Transformer):
    """Test helper that forwards bytes in fixed chunks and records pull count."""

    def __init__(self, chunk_size: int = 2):
        super().__init__()
        self.chunk_size = chunk_size
        self.read_calls = 0

    def _fill_buffer(self) -> bytes:
        self.read_calls += 1
        if self.source is None:
            return b""
        return self.source.read(self.chunk_size)


class _RecordingObserver(Observer):
    """Test helper that records observed chunks while forwarding to its sink."""

    def __init__(self):
        super().__init__()
        self.chunks: list[bytes] = []

    def observe(self, chunk: bytes) -> None:
        self.chunks.append(chunk)


class TestPipeAssociativity:
    def test_read_streams_associative(self):
        x = ReadStream(BytesIO(b"abcdef"))
        y = _CountingChunkTransformer(chunk_size=2)
        z = _CountingChunkTransformer(chunk_size=2)
        left = x | (y | z)

        assert left is z
        assert left.source is y
        assert y.source is x

        x2 = ReadStream(BytesIO(b"abcdef"))
        y2 = _CountingChunkTransformer(chunk_size=2)
        z2 = _CountingChunkTransformer(chunk_size=2)
        right = x2 | y2 | z2

        assert right is z2
        assert right.source is y2
        assert y2.source is x2

        assert left.read() == right.read() == b"abcdef"

        assert y.read_calls > 0
        assert z.read_calls > 0
        assert y.read_calls == y2.read_calls
        assert z.read_calls == z2.read_calls

    def test_write_streams_associative(self):
        x = _RecordingObserver()
        y = _RecordingObserver()
        z = BytesIO()
        left = x | (y | z)
        left.write(b"ab")
        left.write(b"c")

        x2 = _RecordingObserver()
        y2 = _RecordingObserver()
        z2 = BytesIO()
        right = x2 | y2 | z2
        right.write(b"ab")
        right.write(b"c")

        assert left is x
        assert right is x2
        assert x.sink is y
        assert y.sink is z
        assert x2.sink is y2
        assert y2.sink is z2
        assert x.chunks == [b"ab", b"c"]
        assert x2.chunks == [b"ab", b"c"]
        assert y.chunks == [b"ab", b"c"]
        assert y2.chunks == [b"ab", b"c"]
        assert z.getvalue() == b"abc"
        assert z2.getvalue() == b"abc"


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
