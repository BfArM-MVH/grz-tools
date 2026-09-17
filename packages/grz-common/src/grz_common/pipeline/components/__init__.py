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

from grz_common.exceptions import UploadError

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


class DataIntegrityError(PipelineError, UploadError):
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

    ``read()`` waits for the first chunk only, then also takes the chunks already queued, up to
    ``size`` bytes, and returns them joined. A reader in another thread, such as grz_check, then
    needs the GIL once per batch instead of once per chunk, and each batch is copied only once.
    ``readinto()`` goes through ``read()``, and ``readall()`` comes from ``io.RawIOBase``.
    """

    def __init__(self, max_queue_size: int = 128) -> None:
        self.queue: queue.Queue[bytes | None] = queue.Queue(maxsize=max_queue_size)
        self.buffer = b""  # the rest of a chunk that did not fit into the last read()
        self.eof = False

    def readable(self) -> bool:
        return True

    def read(self, size: int | None = -1) -> bytes:
        if size is None or size < 0:
            return self.readall()
        parts: list[bytes] = []
        n = 0
        while n < size:
            if self.buffer:
                chunk = self.buffer
                self.buffer = b""
            elif self.eof:
                break
            else:
                try:
                    item = self.queue.get(block=not parts)
                except queue.Empty:
                    break
                if item is None:
                    self.eof = True
                    break
                chunk = item
            room = size - n
            if len(chunk) > room:
                self.buffer = chunk[room:]  # the rest waits for the next read()
                chunk = chunk[:room]
            parts.append(chunk)
            n += len(chunk)
        return b"".join(parts)

    def readinto(self, buffer: Buffer) -> int:
        view = memoryview(buffer).cast("B")
        data = self.read(len(view))
        view[: len(data)] = data
        return len(data)


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
