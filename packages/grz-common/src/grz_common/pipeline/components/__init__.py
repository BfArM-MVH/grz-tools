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
from typing import Any, Protocol, runtime_checkable

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


class Pipeable:
    """Mixin to enable '|' (chaining) and '>>' (execution) operators for streaming."""

    def __or__(self, other: Any) -> Any:
        """
        Piping operator for chaining components.

        1. self: Pipeable | Class -> Class(self)
        2. self: Readable | other: ReadStream  -> other.set_source(self)
        3. self: WriteStream | other: Writable -> self.set_sink(other)
        """
        if isinstance(other, type):
            return other(self)

        if isinstance(self, Readable) and self.readable() and hasattr(other, "source"):
            other.source = self
            return other

        if isinstance(self, WriteStream) and self.writable() and isinstance(other, Writable):
            self.sink = other
            return self

        raise TypeError(f"Operator '|' expects a Transformer or Observer, got {type(other)}")

    def __rshift__(self, other: Writable) -> Writable:
        """
        Redirection operator for driving the pipeline into a sink.
        Usage: (source | transform) >> sink
        """
        if not isinstance(other, Writable):
            raise TypeError(f"Operator '>>' expects a Writable sink, got {type(other)}")

        if not isinstance(self, Readable):
            raise TypeError(f"Cannot drive pipeline: {type(self)} is not Readable.")

        try:
            shutil.copyfileobj(self, other, length=READ_CHUNK_SIZE)
        finally:
            # Ensure close propagates down the whole chain to join threads and flush buffers
            if hasattr(self, "close"):
                self.close()
            # Sinks like open files or DevNullSink need closing to ensure flush
            other.close()
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
        self._source = source

    def read(self, size: int | None = -1) -> bytes:
        if self._source is None:
            raise StreamConfigurationError("Stream source not set. Use '|' to attach a source.")
        return self._source.read(size)

    def close(self) -> None:
        if not self.closed:
            if self._source:
                with contextlib.suppress(Exception):
                    self._source.close()
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

        if not self._output_buffer:
            chunk = self._fill_buffer()
            if not chunk:
                return b""
            if target_size == -1 or len(chunk) == target_size:
                return chunk
            self._output_buffer.extend(chunk)

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
            self._sink = sink

    def write(self, data: Buffer) -> int:
        if self._sink is None:
            raise StreamConfigurationError("Stream sink not set.")
        return self._sink.write(data)

    def close(self) -> None:
        if not self.closed:
            if self._sink:
                self._sink.close()
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


class Tee(ReadStream):
    """
    Branches the stream to an Observer.
    Supports asynchronous background threads or synchronous execution.
    """

    def __init__(self, observer: Observer, max_queue_size: int = 128, threaded: bool = False):
        super().__init__(None)
        self.observer = observer
        self.max_queue_size = max_queue_size
        self._threaded = threaded
        self._queue: queue.Queue[bytes | None] | None = None
        self._thread: threading.Thread | None = None
        self._exc: Exception | None = None

    @property
    def source(self) -> Readable | None:
        return self._source

    @source.setter
    def source(self, source: Readable) -> None:
        super(Tee, type(self)).source.fset(self, source)  # type: ignore[attr-defined]

        if source is not None and self._threaded:
            self._queue = queue.Queue(maxsize=self.max_queue_size)
            self._thread = threading.Thread(target=self._worker, daemon=True)
            self._thread.start()

    def read(self, size: int | None = -1) -> bytes:
        if self._threaded and self._exc:
            raise self._exc

        if self._source is None:
            raise RuntimeError("Stream source not set")

        chunk = self._source.read(size)

        if self._threaded and self._queue:
            if chunk:
                self._queue.put(chunk, block=True)
            else:
                self._queue.put(None)
        elif chunk:
            self.observer.write(chunk)

        return chunk

    def _worker(self) -> None:
        try:
            if not self._queue:
                return
            while True:
                chunk = self._queue.get()
                if chunk is None:
                    self._queue.task_done()
                    break
                self.observer.write(chunk)
                self._queue.task_done()
        except Exception as e:
            self._exc = e

    def close(self) -> None:
        # close upstream sources first
        super().close()

        # shut down the background worker gracefully
        if self._threaded and self._thread and self._thread.is_alive() and self._queue:
            with contextlib.suppress(queue.Full):
                self._queue.put_nowait(None)
            self._thread.join()

        if hasattr(self.observer, "close"):
            self.observer.close()

        if self._exc:
            raise self._exc


class TqdmObserver(Observer):
    def __init__(self, pbar: Any | list[Any]):
        super().__init__()
        self.pbars = pbar if isinstance(pbar, list) else [pbar]
        self._lock = threading.Lock()

    def observe(self, chunk: bytes) -> None:
        n = len(chunk)
        with self._lock:
            for pbar in self.pbars:
                pbar.update(n)


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
