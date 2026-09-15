"""
Base classes and implementations for pipeline components.
Uses '|' (or) for chaining and '>>' (rshift) for execution.
"""

import abc
import contextlib
import io
import logging
import queue
import shutil
import threading
from collections.abc import Buffer
from types import TracebackType
from typing import Any, Protocol, Self, runtime_checkable

log = logging.getLogger(__name__)

READ_CHUNK_SIZE = 8 * 1024 * 1024


class PipelineError(Exception):
    """Base exception for pipeline errors."""

    def __init__(self, message: str, stage: str | None = None, cause: Exception | None = None):
        self.stage = stage
        self.cause = cause
        msg = f"[{stage}] {message}" if stage else message
        if cause:
            msg += f" (Caused by: {type(cause).__name__}: {cause})"
        super().__init__(msg)


class StreamStateError(PipelineError):
    """Raised when an operation is attempted on a closed or unusable stream."""

    pass


class StreamConfigurationError(PipelineError):
    """Raised when the pipeline setup is invalid (e.g., missing source/sink)."""

    pass


class DataValidationError(PipelineError):
    """Raised when data content fails validation (checksum, FASTQ/BAM format, etc.)."""

    pass


class DataIntegrityError(PipelineError):
    """Raised when transfer integrity fails (e.g., S3 ETag mismatch)."""

    pass


@runtime_checkable
class Readable(Protocol):
    """Protocol for objects that can be read from."""

    def read(self, size: int | None = -1) -> bytes: ...
    def readable(self) -> bool: ...
    @property
    def closed(self) -> bool: ...
    def close(self) -> None: ...


@runtime_checkable
class Writable(Protocol):
    """Protocol for objects that can be written to."""

    def write(self, data: Buffer) -> int: ...
    def writable(self) -> bool: ...
    @property
    def closed(self) -> bool: ...
    def close(self) -> None: ...
    # Sinks are driven as context managers by '>>' so they can finalize or abort.
    # Signature matches io.IOBase so io-based sinks satisfy the protocol cleanly.
    def __enter__(self) -> Self: ...
    def __exit__(
        self, exc_type: type[BaseException] | None, exc_val: BaseException | None, exc_tb: TracebackType | None
    ) -> None: ...


class Pipeable:
    """Mixin to enable '|' (chaining) and '>>' (execution) operators for streaming."""

    def __or__(self, other: Any) -> Any:
        """
        Piping operator for chaining components.

        0. self | None -> self
        1. self | callable -> callable(self), e.g. a class or StreamMetricsRegistry.measure(...)
        2. self: Readable | other: ReadStream  -> other, with other.source = self
            -> other will read from self
        3. self: WriteStream | other: Writable -> self, with self.sink = other
            -> self will write to other
        """
        if other is None:
            return self

        if callable(other):
            return other(self)

        if isinstance(self, Readable) and self.readable() and isinstance(other, ReadStream):
            other.source = self
            return other

        if isinstance(self, WriteStream) and isinstance(other, Writable) and other.writable():
            self.sink = other
            return self

        raise TypeError(f"Operator '|' expects a pipeable object or a callable, got {type(other)}")

    def __rshift__(self, other: Writable) -> Writable:
        """
        Redirection operator for driving the pipeline into a sink.
        Usage: (source | transform) >> sink
        """
        if not isinstance(other, Writable):
            raise TypeError(f"Operator '>>' expects a Writable sink, got {type(other)}")

        if not isinstance(self, Readable):
            raise TypeError(f"Cannot drive pipeline: {type(self)} is not Readable.")

        # Drive everything through the destination's context manager: it finalizes the write
        # on a clean exit and aborts it on error. Close the source inside the block so that a
        # validation failure there also triggers the abort, rather than leaving a half-uploaded
        # object behind. If streaming itself failed, report that error, not a follow-up one
        # from closing the stages.
        with other:
            try:
                shutil.copyfileobj(self, other, length=READ_CHUNK_SIZE)
            except BaseException:
                with contextlib.suppress(Exception):
                    self.close()
                raise
            self.close()
        return other


class ReadStream(io.BufferedIOBase, Pipeable):
    """Wraps a source input."""

    def __init__(self, source: Readable | None = None):
        super().__init__()
        self._source: Readable | None = source

    def readable(self) -> bool:
        return True

    @property
    def source(self) -> Readable | None:
        return self._source

    @source.setter
    def source(self, source: Readable) -> None:
        if self.closed:
            raise StreamStateError("Cannot set source on a closed stream")

        if not isinstance(source, Readable):
            raise TypeError(f"Source must be a readable object. Got: {type(source).__name__}")

        if self._source is None:
            self._source = source
        elif isinstance(self._source, ReadStream):
            self._source.source = source
        else:
            raise TypeError(f"Cannot set source: {type(self._source)} is not a ReadStream")

    def read(self, size: int | None = -1) -> bytes:
        if self._source is None:
            raise StreamConfigurationError("Stream source not set. Use '|' to attach a source.")
        return self._source.read(size)

    def close(self) -> None:
        if self.closed:
            return
        try:
            if self._source:
                self._source.close()
        finally:
            super().close()


class Transformer(ReadStream, metaclass=abc.ABCMeta):
    """
    Reads from upstream, transforms data, yields to downstream.
    """

    def __init__(self, source: Readable | None = None):
        super().__init__(source)
        self._output_buffer = bytearray()

    @abc.abstractmethod
    def _fill_buffer(self) -> bytes:
        """
        Override this: Read from self.source, transform, return bytes.
        Return empty bytes b"" on EOF.
        """
        raise NotImplementedError

    def read(self, size: int | None = -1) -> bytes:
        target_size = size if size is not None else -1

        while target_size == -1 or len(self._output_buffer) < target_size:
            chunk = self._fill_buffer()
            if not chunk:
                break
            self._output_buffer.extend(chunk)

        limit = len(self._output_buffer) if target_size == -1 else min(len(self._output_buffer), target_size)
        ret = self._output_buffer[:limit]
        del self._output_buffer[:limit]
        return bytes(ret)


class WriteStream(io.BufferedIOBase, Pipeable):
    """Wraps a sink output."""

    def __init__(self, sink: Writable | None = None):
        super().__init__()
        self._sink: Writable | None = sink

    def writable(self) -> bool:
        return True

    @property  # type: ignore[override]
    def sink(self) -> Writable | None:
        return self._sink

    @sink.setter
    def sink(self, sink: Writable) -> None:
        if self.closed:
            raise StreamStateError("Cannot set sink on a closed stream")

        if not isinstance(sink, Writable):
            raise TypeError(f"Sink must be a writable object with a 'write' method. Got: {type(sink).__name__}")

        if self._sink is None:
            self._sink = sink
        elif isinstance(self._sink, WriteStream):
            self._sink.sink = sink
        else:
            raise TypeError(f"Cannot set sink: {type(self._sink)} is not a WriteStream")

    def write(self, data: Buffer) -> int:
        if self._sink is None:
            raise StreamConfigurationError("Stream sink not set.")
        return self._sink.write(data)

    def close(self) -> None:
        if self.closed:
            return
        try:
            if self._sink:
                self._sink.close()
        finally:
            super().close()


class Observer(WriteStream, metaclass=abc.ABCMeta):
    """Accepts data via write(), processes it, and pushes to next observer (if any)."""

    def write(self, data: Buffer) -> int:
        # use memoryview to handle the abstract Buffer type
        mv = memoryview(data)
        # observe protocol demands bytes
        self.observe(mv.tobytes())

        if self._sink:
            return self._sink.write(data)

        return len(mv)

    @abc.abstractmethod
    def observe(self, chunk: bytes) -> None:
        raise NotImplementedError()


class Metrics(Protocol):
    @property
    def metrics(self) -> dict[str, Any]: ...


class ObserverWithMetrics(Observer, Metrics, metaclass=abc.ABCMeta):
    pass


class PushToPullAdapter(io.RawIOBase):
    """
    File-like adapter that bridges pipeline push operations to pull ones.

    Producers put byte chunks on ``queue`` and ``None`` at the end of the stream.
    ``read()`` and ``readall()`` come from ``io.RawIOBase`` via ``readinto()``.

    ``readinto()`` waits for the first chunk only, then also takes the chunks already queued, up
    to the size of the buffer. A reader in another thread, such as grz_check, then needs the GIL
    once per batch instead of once per chunk.
    """

    def __init__(self, max_queue_size: int = 128) -> None:
        self.queue: queue.Queue[bytes | None] = queue.Queue(maxsize=max_queue_size)
        self.buffer = bytearray()
        self.eof = False

    def readable(self) -> bool:
        return True

    def readinto(self, buffer: Buffer) -> int:
        view = memoryview(buffer).cast("B")
        while len(self.buffer) < len(view) and not self.eof:
            try:
                chunk = self.queue.get(block=not self.buffer)
            except queue.Empty:
                break
            if chunk is None:
                self.eof = True
            else:
                self.buffer.extend(chunk)

        n = min(len(view), len(self.buffer))
        view[:n] = self.buffer[:n]
        del self.buffer[:n]
        return n


class Tee(ReadStream):
    """Branches the stream to an observer."""

    def __init__(self, observer: Writable):
        super().__init__(None)
        self.observer = observer

    def read(self, size: int | None = -1) -> bytes:
        chunk = super().read(size)
        if chunk:
            self.observer.write(chunk)
        return chunk

    def close(self) -> None:
        if self.closed:
            return
        # close upstream sources first
        try:
            super().close()
        finally:
            self.observer.close()


class TqdmObserver(Observer):
    # One lock for all instances: threads share progress bars (e.g. the total bar of a
    # thread pool), and tqdm's update() does not lock its counter. Other updates of a
    # shared bar must take it too.
    lock = threading.Lock()

    def __init__(self, pbar: Any | list[Any]):
        super().__init__()
        self.pbars = pbar if isinstance(pbar, list) else [pbar]

    def observe(self, chunk: bytes) -> None:
        n = len(chunk)
        with self.lock:
            for pbar in self.pbars:
                pbar.update(n)


class ThreadedObserver(Observer):
    """
    Runs an observer in its own worker thread, so a slow observer does not hold up the stream.

    Chunks travel to the worker through a queue that holds at most ``max_queue_size`` chunks.
    If the observer fails, the next write or ``close()`` raises its error; neither waits on a
    queue that nobody empties anymore. ``close()`` waits until the worker has written every
    chunk, then closes the observer.

    Usage::

        pipeline |= Tee(ThreadedObserver(ChecksumValidator(expected_checksum=...)))
    """

    def __init__(self, observer: Writable, max_queue_size: int = 8) -> None:
        super().__init__()
        self.observer = observer
        self._queue: queue.Queue[bytes | None] = queue.Queue(maxsize=max_queue_size)
        self._thread: threading.Thread | None = None
        self._exc: Exception | None = None

    def observe(self, chunk: bytes) -> None:
        if self.closed:
            raise ValueError("I/O operation on closed file.")
        if self._thread is None:
            # start on the first chunk, not when the pipeline is built
            self._thread = threading.Thread(target=self._work, daemon=True)
            self._thread.start()
        self._put(chunk)

    def _work(self) -> None:
        try:
            while (chunk := self._queue.get()) is not None:
                self.observer.write(chunk)
        except Exception as e:
            self._exc = e

    def _put(self, item: bytes | None) -> None:
        """Hand an item to the worker, without waiting forever once the worker has stopped."""
        while True:
            self._raise_if_failed()
            try:
                self._queue.put(item, timeout=0.1)
                return
            except queue.Full:
                if self._thread is not None and not self._thread.is_alive():
                    self._raise_if_failed()
                    raise PipelineError("Observer thread stopped unexpectedly", stage=self.__class__.__name__) from None

    def _raise_if_failed(self) -> None:
        if self._exc is not None:
            raise self._exc

    def close(self) -> None:
        if self.closed:
            return
        try:
            try:
                if self._thread is not None:
                    self._put(None)  # end marker
                    self._thread.join()
                    self._raise_if_failed()
            except BaseException:
                # report the worker's error, not a follow-up one from closing the observer
                with contextlib.suppress(Exception):
                    self.observer.close()
                raise
            self.observer.close()
        finally:
            super().close()


class DevNullSink(io.BufferedIOBase, Writable):
    """Sink that discards all data."""

    def writable(self) -> bool:
        return True

    def write(self, data: Buffer) -> int:
        return len(memoryview(data))

    @property
    def closed(self) -> bool:
        return super().closed

    def close(self) -> None:
        super().close()
