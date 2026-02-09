"""Base classes and interfaces for pipeline components."""

from __future__ import annotations

import logging
from abc import abstractmethod
from collections.abc import Iterator
from dataclasses import dataclass, field
from typing import Any

log = logging.getLogger(__name__)


class PipelineError(Exception):
    """Base exception for pipeline errors."""

    def __init__(self, message: str, stage: str | None = None, cause: Exception | None = None):
        self.stage = stage
        self.cause = cause
        super().__init__(f"[{stage}] {message}" if stage else message)


@dataclass
class PipelineContext:
    """
    Shared context passed through the pipeline.

    Allows stages to share state and collect results.
    """

    # Accumulated errors from all stages
    errors: list[str] = field(default_factory=list)

    # Metrics and counters
    bytes_read: int = 0
    bytes_written: int = 0
    bytes_decrypted: int = 0

    # Checksums (populated by checksummer stages)
    checksums: dict[str, str] = field(default_factory=dict)

    # Validation state
    validation_passed: bool = True

    # Arbitrary stage-specific data
    metadata: dict[str, Any] = field(default_factory=dict)

    def add_error(self, error: str) -> None:
        """Add an error to the context."""
        self.errors.append(error)
        self.validation_passed = False

    def has_errors(self) -> bool:
        """Check if any errors have been recorded."""
        return len(self.errors) > 0


class StreamStage:
    """
    Base class for all pipeline stages.

    Each stage has:
    - A name for identification and logging
    - Access to the shared pipeline context
    - Lifecycle methods (initialize, finalize)

    Supports context manager protocol for automatic cleanup.
    """

    def __init__(self, name: str | None = None):
        self._name = name or self.__class__.__name__
        self._context: PipelineContext | None = None
        self._log = log.getChild(self._name)
        self._initialized = False

    @property
    def name(self) -> str:
        """Return the stage name."""
        return self._name

    @property
    def context(self) -> PipelineContext:
        """Return the pipeline context."""
        if self._context is None:
            raise RuntimeError(f"Stage {self._name} not initialized with context")
        return self._context

    def initialize(self, context: PipelineContext) -> None:
        """
        Initialize the stage with a pipeline context.

        Override to perform setup operations.
        """
        self._context = context
        self._initialized = True

    def finalize(self) -> None:
        """
        Finalize the stage.

        Override to perform cleanup and collect final results.
        Called after all data has been processed.
        """
        pass

    def abort(self) -> None:
        """
        Abort the stage due to an error.

        Override to perform cleanup on failure.
        """
        pass

    def __enter__(self) -> StreamStage:
        """Enter context manager."""
        return self

    def __exit__(self, exc_type: Any, exc_val: Any, exc_tb: Any) -> None:
        """Exit context manager, aborting on error or finalizing on success."""
        if exc_type is not None:
            self.abort()
        elif self._initialized:
            self.finalize()


class StreamSource(StreamStage):
    """
    A source that produces data for the pipeline.

    Examples: S3 downloader, file reader
    """

    @abstractmethod
    def read(self, size: int = -1) -> bytes:
        """
        Read data from the source.

        :param size: Number of bytes to read (-1 for all available)
        :returns: Bytes read, empty bytes when exhausted
        """
        pass

    def __iter__(self) -> Iterator[bytes]:
        """Iterate over chunks from the source."""
        while chunk := self.read(65536):  # Default chunk size
            yield chunk

    def iter_chunks(self, chunk_size: int = 65536) -> Iterator[bytes]:
        """Iterate over chunks of specified size."""
        while chunk := self.read(chunk_size):
            yield chunk

    @property
    def content_length(self) -> int | None:
        """Return content length if known, None otherwise."""
        return None


class StreamTransformer(StreamStage):
    """
    A stage that transforms data passing through.

    Examples: Decryptor, Decompressor, Encryptor

    Transformers may buffer data internally and produce output
    asynchronously (e.g., encryption produces output in segment-sized chunks).
    """

    @abstractmethod
    def process(self, data: bytes) -> bytes:
        """
        Process a chunk of input data.

        :param data: Input bytes to transform
        :returns: Transformed output bytes (may be empty if buffering)
        """
        pass

    def flush(self) -> bytes:
        """
        Flush any buffered data.

        Called during finalization to ensure all data is processed.

        :returns: Any remaining output data
        """
        return b""


class StreamObserver(StreamStage):
    """
    A stage that observes data without modification.

    Examples: Checksummer, Line counter, Validator

    Observers receive data and update their internal state,
    but always return the data unchanged.
    """

    @abstractmethod
    def observe(self, data: bytes) -> None:
        """
        Observe a chunk of data.

        :param data: Bytes to observe (will not be modified)
        """
        pass

    def get_result(self) -> Any:
        """
        Get the observation result.

        Override to return accumulated results (e.g., checksum, count).
        """
        return None


class StreamSink(StreamStage):
    """
    A sink that consumes data from the pipeline.

    Examples: S3 uploader, file writer
    """

    @abstractmethod
    def write(self, data: bytes) -> int:
        """
        Write data to the sink.

        :param data: Bytes to write
        :returns: Number of bytes written
        """
        pass

    @property
    def bytes_written(self) -> int:
        """Return total bytes written."""
        return 0
