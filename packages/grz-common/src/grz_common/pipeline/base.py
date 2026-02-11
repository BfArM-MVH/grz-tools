"""Base classes and interfaces for pipeline components."""

from __future__ import annotations

import io
import logging
from abc import ABCMeta, abstractmethod

log = logging.getLogger(__name__)


class PipelineError(Exception):
    """Base exception for pipeline errors."""

    def __init__(self, message: str, stage: str | None = None, cause: Exception | None = None):
        self.stage = stage
        self.cause = cause
        super().__init__(f"[{stage}] {message}" if stage else message)


class StreamWrapper(io.BufferedIOBase):
    """
    Base class that simply wraps another stream.
    Delegates close/flush/context to the source.
    """

    def __init__(self, source: io.BufferedIOBase):
        self.source = source

    def readable(self) -> bool:
        return True

    def close(self):
        # We assume ownership of the source stream and close it too.
        if not self.closed:
            self.source.close()
            super().close()


class TransformStream(StreamWrapper, metaclass=ABCMeta):
    """
    Base class for streams that MODIFY data (Encryption, Decryption).
    Guarantees full reads (up to 'size') unless EOF is reached.
    """

    def __init__(self, source: io.BufferedIOBase):
        super().__init__(source)
        self._output_buffer = bytearray()
        self._eof = False

    @abstractmethod
    def _fill_buffer(self) -> None:
        """
        Logic to read from source, transform, and append to self._output_buffer.
        Must set self._eof = True when source is exhausted.
        """
        raise NotImplementedError

    def read(self, size: int = -1) -> bytes:
        # 1. Handle "Read All" (size=-1) or "Read until EOF"
        if size == -1:
            while not self._eof:
                self._fill_buffer()
            ret = bytes(self._output_buffer)
            self._output_buffer.clear()
            return ret

        # 2. Serve from buffer if we already have enough
        if len(self._output_buffer) >= size:
            ret = self._output_buffer[:size]
            self._output_buffer = self._output_buffer[size:]
            return bytes(ret)

        # 3. Fill buffer until we have enough data or hit EOF
        while len(self._output_buffer) < size and not self._eof:
            self._fill_buffer()

        # 4. Return what we have (may be short if EOF hit)
        limit = min(len(self._output_buffer), size)
        ret = self._output_buffer[:limit]
        self._output_buffer = self._output_buffer[limit:]
        return bytes(ret)


class ObserverStream(StreamWrapper, metaclass=ABCMeta):
    """
    Base class for streams that INSPECT data (Validation, Hashing).
    Pass-through: Does NOT modify data.
    """

    @abstractmethod
    def observe(self, chunk: bytes) -> None:
        raise NotImplementedError

    def read(self, size: int = -1) -> bytes:
        chunk = self.source.read(size)
        if chunk:
            self.observe(chunk)
        return chunk
