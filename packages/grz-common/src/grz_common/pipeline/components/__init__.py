"""Base classes and implementations for pipeline components."""

import contextlib
import io
import threading
from abc import ABCMeta, abstractmethod
from typing import Any

READ_CHUNK_SIZE = 1 * 1024 * 1024


class PipelineError(Exception):
    """Base exception for pipeline errors."""

    def __init__(self, message: str, stage: str | None = None, cause: Exception | None = None):
        self.stage = stage
        self.cause = cause
        super().__init__(f"[{stage}] {message}" if stage else message)


class StreamWrapper(io.BufferedIOBase):
    def __init__(self, source: io.BufferedIOBase):
        self.source = source

    def readable(self) -> bool:
        return True

    def close(self):
        if not self.closed:
            with contextlib.suppress(Exception):
                self.source.close()
            super().close()


class ObserverStream(StreamWrapper, metaclass=ABCMeta):
    @abstractmethod
    def observe(self, chunk: bytes) -> None:
        raise NotImplementedError

    def read(self, size: int | None = -1) -> bytes:
        chunk = self.source.read(size)
        if chunk:
            self.observe(chunk)
        return chunk


class TqdmObserver(ObserverStream):
    def __init__(self, source: io.BufferedIOBase, pbar: Any | list[Any]):
        super().__init__(source)
        self.pbars = pbar if isinstance(pbar, list) else [pbar]
        self._lock = threading.Lock()

    def observe(self, chunk: bytes) -> None:
        n = len(chunk)
        with self._lock:
            for pbar in self.pbars:
                pbar.update(n)


class TransformStream(StreamWrapper, metaclass=ABCMeta):
    """
    Base class for streams that MODIFY data.
    """

    def __init__(self, source: io.BufferedIOBase):
        super().__init__(source)
        self._output_buffer = bytearray()
        self._eof = False

    @abstractmethod
    def _fill_buffer(self) -> bytes:
        """
        Read and transform data from the source stream.

        Subclasses must implement this method to read from the source, apply any transformation, and return the processed bytes.

        Example for a simple passthrough implementation::

            def _fill_buffer(self) -> bytes:
                return self.source.read(READ_CHUNK_SIZE)

        :returns: Transformed data chunk, or empty bytes when EOF is reached.
        """
        raise NotImplementedError

    def read(self, size: int | None = -1) -> bytes:
        target_size = size if size is not None else -1
        read_all = target_size == -1

        # Keep filling the output buffer until we have enough data to return or reach EOF
        while (read_all or len(self._output_buffer) < target_size) and (chunk := self._fill_buffer()):
            self._output_buffer.extend(chunk)

        # Return up to target_size bytes from the output buffer
        limit = min(len(self._output_buffer), target_size)
        ret = self._output_buffer[:limit]
        del self._output_buffer[:limit]
        return bytes(ret)


class ValidatingStream(ObserverStream, metaclass=ABCMeta):
    """
    Base class for streams that perform validation.
    Automatically runs validate() when the context manager exits.
    """

    @abstractmethod
    def validate(self) -> None:
        """
        Run final validity checks.
        Should raise an exception if validation fails.
        """
        pass

    @property
    def metrics(self) -> dict[str, Any]:
        return {}

    def close(self):
        """
        Trigger validation on close if it hasn't happened yet.
        """
        if not self.closed:
            try:
                self.validate()
            except Exception:
                # if validation fails, we still want to try closing the underlying stream,
                # but we must re-raise the validation error.
                super().close()
                raise
            super().close()

    def __exit__(self, exc_type, exc_val, exc_tb):
        # if an exception already occurred in the block, we don't want to mask it
        # with a validation error unless we want to enforce validation even on failure.
        if exc_type is None:
            self.validate()
        super().__exit__(exc_type, exc_val, exc_tb)
