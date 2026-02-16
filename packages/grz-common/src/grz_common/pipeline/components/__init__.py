"""
Base classes and implementations for pipeline components.
Uses '|' (or) for chaining and '>>' (rshift) for execution.
"""

import contextlib
import io
import logging
import queue
import shutil
import threading
from io import IOBase
from typing import Any, Protocol, runtime_checkable

log = logging.getLogger(__name__)

READ_CHUNK_SIZE = 8 * 1024 * 1024


class PipelineError(Exception):
    """Base exception for pipeline errors."""

    def __init__(self, message: str, stage: str | None = None, cause: Exception | None = None):
        self.stage = stage
        self.cause = cause
        super().__init__(f"[{stage}] {message}" if stage else message)


@runtime_checkable
class Writable(Protocol):
    def write(self, data: bytes) -> int: ...
    def close(self) -> None: ...


class Pipeable:
    """Mixin to enable '|' (chaining) and '>>' (execution) operators for streaming."""

    def __or__(self, other):
        """
        Piping operator for chaining components.

        1. stream | TransformerClass -> TransformerClass(stream)
        2. stream | transformer_instance -> transformer_instance.source = stream
        3. observer | observer -> Links push-based observers in a side-channel.
        """
        if isinstance(other, type):
            return other(self)

        if isinstance(self, Observer) and isinstance(other, Observer):
            self.set_sink(other)
            return self

        if isinstance(other, (Transformer, Tee)):
            other.source = self
            return other

        raise TypeError(f"Operator '|' expects a Transformer or Observer, got {type(other)}")

    def __rshift__(self, other: Writable):
        """
        Redirection operator for driving the pipeline into a sink.
        Usage: (source | transform) >> sink
        """
        if not hasattr(other, "write") or not callable(other.write):
            raise TypeError(f"Operator '>>' expects a Writable sink, got {type(other)}")

        if not hasattr(self, "read"):
            raise TypeError(f"Cannot drive pipeline: {type(self)} is not readable.")

        try:
            # Drive the pull-based pipeline into the sink
            shutil.copyfileobj(self, other, length=READ_CHUNK_SIZE)
        finally:
            # Ensure close propagates down the whole chain to join threads and flush buffers
            if hasattr(self, "close"):
                self.close()
            if hasattr(other, "close"):
                # Sinks like open files or DevNullSink need closing to ensure flush
                other.close()
        return other


class Stream(io.BufferedIOBase, Pipeable):
    """Wraps a source input."""

    def __init__(self, source: io.BufferedIOBase | None = None):
        super().__init__()
        self._source = source

    @property
    def source(self):
        return self._source

    @source.setter
    def source(self, value):
        self._source = value

    def read(self, size: int | None = -1) -> bytes:
        if self.source is None:
            raise RuntimeError("Stream source not set. Use '|' to attach a source.")
        return self.source.read(size)

    def close(self):
        if not self.closed:
            if self._source:
                with contextlib.suppress(Exception):
                    self._source.close()
            super().close()


class Transformer(Stream):
    """
    Reads from upstream, transforms data, yields to downstream.
    """

    def __init__(self, source: io.BufferedIOBase | None = None):
        super().__init__(source)
        self._output_buffer = bytearray()

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


class Observer(Pipeable):
    """
    Accepts data via write(), processes it, and pushes to next observer (if any).
    """

    def __init__(self):
        self.sink: Writable | None = None

    def set_sink(self, sink: Writable):
        self.sink = sink
        return self

    def write(self, data: bytes) -> int:
        self.observe(data)
        if self.sink:
            self.sink.write(data)
        return len(data)

    def observe(self, chunk: bytes) -> None:
        """Subclasses implement logic here."""
        pass

    def close(self):
        if self.sink:
            self.sink.close()


class Metrics(Protocol):
    @property
    def metrics(self) -> dict[str, Any]: ...


class ObserverWithMetrics(Observer, Metrics):
    pass


class Tee(Stream):
    """
    Branches the stream to an Observer.
    Supports asynchronous background threads or synchronous execution.
    """

    def __init__(self, observer: Observer, max_queue_size: int = 128, threaded: bool = False):
        self.observer = observer
        self.max_queue_size = max_queue_size
        self._threaded = threaded
        self._queue: queue.Queue | None = None
        self._thread: threading.Thread | None = None
        self._exc: Exception | None = None
        super().__init__(None)

    @Stream.source.setter  # type: ignore[attr-defined]
    def source(self, value):
        self._source = value
        if value is not None and self._threaded:
            self._queue = queue.Queue(maxsize=self.max_queue_size)
            self._thread = threading.Thread(target=self._worker, daemon=True)
            self._thread.start()

    def read(self, size: int | None = -1) -> bytes:
        chunk = self._source.read(size)  # type: ignore[union-attr]
        if chunk:
            if self._exc:
                raise self._exc

            if self._threaded:
                self._queue.put(chunk, block=True)  # type: ignore[union-attr]
            else:
                try:
                    self.observer.write(chunk)
                except Exception as e:
                    self._exc = e
                    raise
        return chunk

    def _worker(self):
        try:
            while True:
                chunk = self._queue.get()  # type: ignore[union-attr]
                if chunk is None:
                    break
                self.observer.write(chunk)
                self._queue.task_done()  # type: ignore[union-attr]
            self.observer.close()
        except Exception as e:
            self._exc = e

    def close(self):
        if self._threaded:
            if self._thread and self._thread.is_alive() and self._queue:
                self._queue.put(None)
                self._thread.join()
        else:
            self.observer.close()

        if self._exc:
            raise self._exc
        super().close()


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


class DevNullSink(Writable, io.BufferedIOBase):
    """Sink that discards all data."""

    def writable(self) -> bool:
        return True

    def write(self, data) -> int:
        return len(data)
