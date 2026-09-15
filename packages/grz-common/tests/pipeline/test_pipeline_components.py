"""Tests for the modular pipeline components."""

import array
import gzip
import hashlib
import io
import threading
import time
from collections.abc import Callable
from io import BytesIO

import pytest
from grz_common.pipeline.components import (
    DataValidationError,
    Observer,
    PushToPullAdapter,
    ReadStream,
    Tee,
    ThreadedObserver,
    Transformer,
)
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


class _SlowObserver(Observer):
    """Test helper that takes a while per chunk and records the chunks and the threads it ran on."""

    def __init__(self, delay: float = 0.0):
        super().__init__()
        self.delay = delay
        self.chunks: list[bytes] = []
        self.threads: set[int] = set()

    def observe(self, chunk: bytes) -> None:
        time.sleep(self.delay)
        self.threads.add(threading.get_ident())
        self.chunks.append(chunk)


class _FailingObserver(Observer):
    """Test helper that fails on its first chunk, after a delay."""

    def __init__(self, delay: float = 0.0):
        super().__init__()
        self.delay = delay

    def observe(self, chunk: bytes) -> None:
        time.sleep(self.delay)
        raise RuntimeError("observer failed")


def _finishes_within(seconds: float, fn: Callable[[], object]) -> None:
    """Run *fn* in a daemon thread; fail the test instead of hanging if it does not return in time."""
    errors: list[BaseException] = []

    def run() -> None:
        try:
            fn()
        except BaseException as e:  # re-raised in the test thread below
            errors.append(e)

    thread = threading.Thread(target=run, daemon=True)
    thread.start()
    thread.join(seconds)
    if thread.is_alive():
        pytest.fail(f"did not finish within {seconds} s")
    if errors:
        raise errors[0]


class TestThreadedObserver:
    """ThreadedObserver runs an observer in a worker thread and never hangs on a failed or slow worker."""

    def test_passes_every_chunk_in_order_on_another_thread(self):
        observer = _SlowObserver()
        pipeline = ReadStream(BytesIO(b"0123456789")) | Tee(ThreadedObserver(observer))

        _finishes_within(5, lambda: pipeline >> BytesIO())

        assert b"".join(observer.chunks) == b"0123456789"
        assert observer.threads
        assert threading.get_ident() not in observer.threads
        assert observer.closed

    def test_failing_observer_does_not_block_the_producer(self):
        """Hang path 1: the worker dies while the queue is full, so nothing empties it anymore."""
        threaded = ThreadedObserver(_FailingObserver(delay=0.2), max_queue_size=1)

        def write_many() -> None:
            for _ in range(100):
                threaded.write(b"x")

        with pytest.raises(RuntimeError, match="observer failed"):
            _finishes_within(5, write_many)

    def test_close_with_a_full_queue_delivers_every_chunk(self):
        """Hang path 2: close() must get its end marker through a full queue."""
        observer = _SlowObserver(delay=0.01)
        threaded = ThreadedObserver(observer, max_queue_size=1)
        for i in range(20):
            threaded.write(bytes([i]))

        _finishes_within(5, threaded.close)

        assert observer.chunks == [bytes([i]) for i in range(20)]
        assert observer.closed

    def test_worker_error_after_the_last_chunk_is_reported_at_close(self):
        threaded = ThreadedObserver(_FailingObserver(delay=0.1))
        threaded.write(b"last chunk")

        with pytest.raises(RuntimeError, match="observer failed"):
            _finishes_within(5, threaded.close)

    def test_observer_error_at_close_fails_the_pipeline(self):
        """A validator that only fails at close, like a checksum mismatch, still fails the pipeline."""
        pipeline = ReadStream(BytesIO(b"data")) | Tee(ThreadedObserver(ChecksumValidator(expected_checksum="0" * 64)))

        with pytest.raises(DataValidationError, match="Checksum mismatch"):
            _finishes_within(5, lambda: pipeline >> BytesIO())
