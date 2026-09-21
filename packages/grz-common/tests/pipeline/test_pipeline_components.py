"""Tests for the modular pipeline components."""

import array
import contextlib
import gc
import gzip
import hashlib
import io
import os
from io import BytesIO

import grz_check
import pytest
from botocore.exceptions import ClientError
from grz_common.exceptions import MissingObjectError
from grz_common.pipeline.components import (
    DataValidationError,
    Observer,
    PushToPullAdapter,
    ReadStream,
    Tee,
    Transformer,
)
from grz_common.pipeline.components.s3 import S3Downloader
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


class _FailingSource(io.RawIOBase):
    """Test helper that returns one chunk, then fails on the next read."""

    def __init__(self):
        self.reads = 0

    def readable(self) -> bool:
        return True

    def read(self, size: int = -1) -> bytes:
        self.reads += 1
        if self.reads > 1:
            raise OSError("source failed")
        return b"some data"


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

        with pytest.raises(DataValidationError, match=r"invalid name prefix|Failed to parse record"):
            with (
                BytesIO(gzip.compress(fastq_content)) as f,
                ReadStream(f) as source,
                FastqValidator(mean_read_length_threshold=12) as validator,
            ):
                source >> validator

    def test_invalid_fastq_seq_qual_length_mismatch(self):
        """Test that sequence and quality length match."""
        fastq_content = b"@read1\n" + b"A" * 16 + b"\n+\n" + b"I" * 15 + b"\n"

        with pytest.raises(
            DataValidationError, match=r"sequence and quality lengths don't match|Failed to parse record"
        ):
            with (
                BytesIO(gzip.compress(fastq_content)) as f,
                ReadStream(f) as source,
                FastqValidator(mean_read_length_threshold=12) as validator,
            ):
                source >> validator

    def test_invalid_fastq_error_reaches_the_writer(self):
        """A FASTQ that fails validation while it is still being written raises the parse error from ``write``."""
        data = gzip.compress(b"not a fastq record\n" + os.urandom(1024 * 1024))
        validator = FastqValidator()
        try:
            # more writes than the validator queues, so the writer is blocked when grz_check stops
            with pytest.raises(DataValidationError, match="invalid name prefix"):
                for start in range(0, len(data), 4096):
                    validator.write(data[start : start + 4096])
        finally:
            with contextlib.suppress(DataValidationError):
                validator.close()


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

        with pytest.raises(DataValidationError, match=r"Checksum mismatch"):
            with (
                BytesIO(data) as f,
                ReadStream(f) as source,
                ChecksumValidator(expected_checksum="0" * 64) as validator,
            ):
                source >> validator


class TestValidatorWorkerDeath:
    """A validator whose grz_check call dies must not let the data pass as validated."""

    def test_a_panicking_grz_check_reaches_the_caller(self, monkeypatch):
        """grz_check is a pyo3 extension, so a Rust panic arrives as a BaseException."""

        class _Panic(BaseException):
            """Stands in for pyo3's PanicException, which derives from BaseException."""

        def panic(*args, **kwargs):
            raise _Panic("simulated panic in grz_check")

        monkeypatch.setattr(grz_check, "validate_fastq", panic)
        # small enough that the validator queue never fills, so nothing else notices the dead worker
        fastq = gzip.compress(b"@read1\nACGT\n+\nIIII\n")

        pipeline = ReadStream(BytesIO(fastq)) | Tee(FastqValidator())

        with pytest.raises(_Panic):
            pipeline >> BytesIO()


class TestCloseErrors:
    """Errors raised in close() by any stage must reach the caller of '>>'."""

    def test_checksum_mismatch_before_later_stages(self):
        """Same shape as the grzctl process chain: validation Tee, transformer, progress Tee."""
        pipeline = (
            ReadStream(BytesIO(b"Test data for checksum validation"))
            | Tee(ChecksumValidator(expected_checksum="0" * 64))
            | _CountingChunkTransformer(chunk_size=4)
            | Tee(_RecordingObserver())
        )

        with pytest.raises(DataValidationError, match=r"Checksum mismatch"):
            pipeline >> BytesIO()

    def test_stream_error_wins_over_close_error(self):
        """When streaming fails, report that failure, not the follow-up failure from close()."""
        pipeline = (
            ReadStream(_FailingSource())
            | Tee(ChecksumValidator(expected_checksum="0" * 64))
            | _CountingChunkTransformer(chunk_size=4)
        )

        with pytest.raises(OSError, match="source failed"):
            pipeline >> BytesIO()


class TestFailedConstruction:
    """A stage that cannot be built must leave nothing behind for the collector to trip over."""

    @pytest.mark.filterwarnings("error::pytest.PytestUnraisableExceptionWarning")
    def test_a_download_that_cannot_open_is_finalized_quietly(self):
        """S3Downloader runs its base constructor first, so the half-built stage still closes."""

        class _MissingObject:
            def get_object(self, **_kwargs):
                raise ClientError({"Error": {"Code": "NoSuchKey"}}, "GetObject")

        with pytest.raises(MissingObjectError):
            S3Downloader(_MissingObject(), "bucket", "key")

        gc.collect()


class TestPushToPullAdapter:
    """The adapter hands queued chunks to a reader, as a raw file object."""

    def test_readinto_fills_any_writable_buffer_bytewise(self):
        """readinto() must count bytes, also for buffers whose items are wider than one byte."""
        adapter = PushToPullAdapter()
        adapter.queue.put(b"abcdefgh")
        adapter.queue.put(None)
        buffer = array.array("i", [0, 0])

        assert adapter.readinto(buffer) == 8
        assert buffer.tobytes() == b"abcdefgh"

    def test_read_returns_queued_chunks_until_the_end_marker(self):
        adapter = PushToPullAdapter()
        for chunk in (b"abc", b"defg", None):
            adapter.queue.put(chunk)

        assert adapter.read(2) == b"ab"
        assert adapter.readall() == b"cdefg"
        assert adapter.read(1) == b""

    def test_read_takes_the_queued_chunks_up_to_the_size_without_waiting(self):
        """read() waits for one chunk only, then takes the chunks already queued."""
        adapter = PushToPullAdapter()
        for chunk in (b"abc", b"defg", b"hij"):
            adapter.queue.put(chunk)

        assert adapter.read(8) == b"abcdefgh"
        # returns the rest at once, although the queue is empty and no end marker came yet
        assert adapter.read(8) == b"ij"
