"""Utility functions and helpers for pipeline operations."""

from __future__ import annotations

import contextlib
import logging
import os
import queue
import signal
import threading
from collections.abc import Callable
from typing import Any

from grz_common.pipeline import StreamStage

log = logging.getLogger(__name__)


class SignalManager:
    """
    Context manager for graceful signal handling.

    Handles SIGINT (Ctrl+C) with two-stage behavior:
    - First signal: Set interrupt flag for graceful shutdown
    - Second signal: Force immediate exit
    """

    def __init__(self, interrupt_event: threading.Event, logger: logging.Logger | None = None):
        """
        Initialize the signal manager.

        :param interrupt_event: Event to set on first interrupt
        :param logger: Optional logger for interrupt messages
        """
        self._interrupt_event = interrupt_event
        self._log = logger or log
        self._original_handler: Any = None

    def __enter__(self) -> SignalManager:
        """Install signal handler."""

        def signal_handler(signum: int, frame: Any) -> None:
            if self._interrupt_event.is_set():
                # second interrupt: force exit
                self._log.warning("Received second interrupt - forcing immediate exit...")

                os._exit(130)  # 128 + SIGINT(2)
            self._log.warning("\nInterrupt received! Gracefully shutting down...")
            self._log.info("Waiting for in-progress operations. Press Ctrl+C again to force exit.")
            self._interrupt_event.set()

        self._original_handler = signal.signal(signal.SIGINT, signal_handler)
        return self

    def __exit__(self, exc_type: Any, exc_val: Any, exc_tb: Any) -> None:
        """Restore original signal handler."""
        if self._original_handler is not None:
            signal.signal(signal.SIGINT, self._original_handler)
            self._original_handler = None


def drain_queue(q: Any) -> None:
    """
    Drain a queue to unblock waiting put() calls.

    :param q: Queue object with get_nowait() method
    """
    try:
        while True:
            q.get_nowait()
    except queue.Empty:
        pass


def safe_join_thread(thread: threading.Thread, timeout: float, logger: logging.Logger | None = None) -> None:
    """
    Safely join a thread with timeout and warning.

    :param thread: Thread to join
    :param timeout: Timeout in seconds
    :param logger: Optional logger for warnings
    """
    _log = logger or log
    thread.join(timeout=timeout)
    if thread.is_alive():
        _log.warning(f"Thread {thread.name} did not finish within {timeout}s")


def abort_all_stages(stages: list[StreamStage], logger: logging.Logger | None = None) -> None:
    """
    Abort all pipeline stages, suppressing exceptions.

    :param stages: List of pipeline stages with abort() method
    :param logger: Optional logger for errors
    """
    _log = logger or log
    for stage in stages:
        with contextlib.suppress(Exception):
            try:
                stage.abort()
            except Exception as e:
                _log.debug(f"Error aborting {stage.name}: {e}")


def finalize_stages_in_order(stages: list[StreamStage], context: Any, logger: logging.Logger | None = None) -> bool:
    """
    Finalize stages in the given order, checking for errors after each.

    This ensures proper cleanup order (e.g., decompressor before validator)
    and allows early abort if validation fails.

    :param stages: List of stages to finalize in order
    :param context: Pipeline context to check for errors
    :param logger: Optional logger
    :returns: True if all stages finalized without errors, False otherwise
    """
    _log = logger or log
    for stage in stages:
        try:
            stage.finalize()
            if hasattr(context, "has_errors") and context.has_errors():
                _log.debug(f"Errors detected after finalizing {stage.name}")
                return False
        except Exception as e:
            _log.error(f"Error finalizing {stage.name}: {e}")
            return False
    return True


def create_worker_thread(target: Callable[[], None], name: str) -> threading.Thread:
    """
    Create and start a named worker thread.

    :param target: Worker function to run
    :param name: Thread name for debugging
    :returns: Started thread
    """
    thread = threading.Thread(target=target, name=name)
    thread.start()
    return thread
