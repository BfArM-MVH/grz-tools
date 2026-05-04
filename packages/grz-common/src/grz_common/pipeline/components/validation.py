import abc
import logging
import queue
import threading
from typing import Any

import grz_check

from . import DataValidationError, ObserverWithMetrics, PipelineError, PushToPullAdapter

log = logging.getLogger(__name__)


class GrzCheckValidator(ObserverWithMetrics, metaclass=abc.ABCMeta):
    """Base class for validators delegating format validation to grz_check."""

    def __init__(self):
        super().__init__()
        self.adapter = PushToPullAdapter()
        self.report: grz_check.ValidationReport = None
        self.exception = None
        self.validation_thread = threading.Thread(target=self._run_validation_thread, daemon=True)
        self.validation_thread.start()

    @abc.abstractmethod
    def _invoke_grz_check(self) -> Any:
        """Invoke the specific grz_check validation function. Returns a ValidationReport."""
        pass

    @abc.abstractmethod
    def _format_error_prefix(self) -> str:
        """Prefix for error messages."""
        pass

    def _run_validation_thread(self):
        try:
            self.report = self._invoke_grz_check()
        except Exception as exception:
            self.exception = exception

    def _enqueue_chunk(self, chunk: bytes | None):
        """Push a chunk to the adapter queue, monitoring thread health."""
        while True:
            if self.exception:
                raise self.exception
            try:
                self.adapter.queue.put(chunk, timeout=0.1)
                break
            except queue.Full as e:
                if not self.validation_thread.is_alive():
                    if self.report and not self.report.is_valid:
                        errors = "; ".join(self.report.errors)
                        prefix = self._format_error_prefix()
                        raise DataValidationError(f"{prefix}: {errors}", stage=self.__class__.__name__) from e
                    raise PipelineError("Validation thread stopped unexpectedly", stage=self.__class__.__name__) from e

    def close(self):
        if self.closed:
            return

        try:
            if not self.exception and self.validation_thread.is_alive():
                self._enqueue_chunk(None)  # signal EOF

            self.validation_thread.join()

            if self.exception:
                if isinstance(self.exception, PipelineError):
                    raise self.exception
                raise DataValidationError(str(self.exception), stage=self.__class__.__name__, cause=self.exception)

            if self.report:
                if not self.report.is_valid:
                    errors = "; ".join(self.report.errors)
                    prefix = self._format_error_prefix()
                    raise DataValidationError(f"{prefix}: {errors}", stage=self.__class__.__name__)
                for warning in self.report.warnings:
                    prefix = self._format_error_prefix().split()[0]
                    log.warning(f"{prefix} Warning: {warning}")

        finally:
            super().close()


class FastqValidator(GrzCheckValidator):
    """Validates a FASTQ stream using grz_check bindings (handles decompression natively)."""

    def __init__(self, mean_read_length_threshold: float | None = None):
        self._threshold = mean_read_length_threshold
        super().__init__()

    def _invoke_grz_check(self) -> Any:
        threshold = int(self._threshold) if self._threshold is not None else None
        return grz_check.validate_fastq(self.adapter, min_mean_read_length=threshold)

    def _format_error_prefix(self) -> str:
        return "FASTQ Invalid"

    def observe(self, chunk: bytes) -> None:
        self._enqueue_chunk(chunk)

    @property
    def metrics(self) -> dict[str, Any]:
        if not self.report:
            return {
                "read_count": 0,
                "mean_read_length": 0.0,
                "total_bases": 0,
                "line_count": 0,
            }

        read_count = self.report.num_records or 0
        mean_len = self.report.mean_read_length or 0.0

        return {
            "read_count": read_count,
            "mean_read_length": mean_len,
            "total_bases": int(read_count * mean_len),
            "line_count": read_count * 4,
        }


class BamValidator(GrzCheckValidator):
    """Validates a BAM stream using grz_check bindings."""

    def __init__(self):
        self._bytes_seen = 0
        super().__init__()

    def _invoke_grz_check(self) -> Any:
        return grz_check.validate_bam(self.adapter)

    def _format_error_prefix(self) -> str:
        return "BAM Invalid"

    def observe(self, chunk: bytes) -> None:
        self._bytes_seen += len(chunk)
        self._enqueue_chunk(chunk)

    @property
    def metrics(self) -> dict[str, Any]:
        return {"size": self._bytes_seen}


class ChecksumValidator(GrzCheckValidator):
    """Validates checksum using grz_check bindings (native speed)."""

    def __init__(self, expected_checksum: str | None = None):
        self._expected = expected_checksum
        self._bytes_seen = 0
        super().__init__()

    def _invoke_grz_check(self) -> Any:
        return grz_check.validate_raw(self.adapter)

    def _format_error_prefix(self) -> str:
        return "Checksum Invalid"

    def observe(self, chunk: bytes) -> None:
        self._bytes_seen += len(chunk)
        self._enqueue_chunk(chunk)

    def close(self):
        super().close()
        if self.report and self._expected and self.report.sha256 != self._expected:
            raise DataValidationError(
                f"Checksum mismatch! Exp: {self._expected}, Got: {self.report.sha256}", stage=self.__class__.__name__
            )

    @property
    def metrics(self) -> dict[str, Any]:
        return {"checksum": self.report.sha256 if self.report else None, "size": self._bytes_seen}
