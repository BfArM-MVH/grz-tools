"""Stop the files that a pool processes in parallel, one level further on each interrupt."""

import logging
import threading
from collections.abc import Iterable
from concurrent.futures import Future, ThreadPoolExecutor, as_completed, wait

log = logging.getLogger(__name__)


class TerminateInterrupt(KeyboardInterrupt):
    """The interrupt that grzctl raises for SIGTERM.

    A supervisor follows SIGTERM with SIGKILL after a grace period, so :func:`wait_for_files`
    handles it like a second Ctrl-C and stops the running files at once.
    """


def wait_for_files(pool: ThreadPoolExecutor, futures: Iterable[Future], stop: threading.Event) -> None:
    """Wait for the files that ``pool`` processes, and stop them one level further on each interrupt.

    1. The first ``KeyboardInterrupt`` starts no queued file, and waits for the running files to finish.
    2. A second one sets ``stop``, so that the running files stop at their next chunk, and waits for that.
    3. A third one stops the waiting and propagates.
       The running files still stop at their next chunk, and the process exits only after that.

    A :class:`TerminateInterrupt` starts at the second level.
    After the first or second level, the first interrupt is raised again.
    Any other error of a file starts no queued file, waits for the running files to finish, and is raised again.
    An interrupt during that wait stops the running files, as the second level does.

    :param pool: The pool that processes the files.
    :param futures: The futures of the files in ``pool``.
    :param stop: The event that the files check before each chunk.
    :raises KeyboardInterrupt: The first interrupt, or the third one.
    """
    futures = list(futures)
    try:
        for future in as_completed(futures):
            future.result()
    except KeyboardInterrupt as interrupt:
        _wait_for_running_files(pool, futures, stop, interrupt)
        raise
    except Exception:
        _wait_for_running_files(pool, futures, stop, None)
        raise


def _wait_for_running_files(
    pool: ThreadPoolExecutor, futures: list[Future], stop: threading.Event, interrupt: KeyboardInterrupt | None
) -> None:
    """Start no queued file, and wait for the running files. Stop them first on a SIGTERM, or on the next interrupt.

    :param interrupt: The interrupt that ended the wait for the files, or ``None`` if a file failed.
    :raises KeyboardInterrupt: An interrupt while the running files stop.
    """
    stopping = isinstance(interrupt, TerminateInterrupt)
    if not stopping:
        try:
            if interrupt is not None:
                log.warning(
                    "Interrupted: no queued file starts, and the running files finish. Press Ctrl-C again to stop them."
                )
            _shut_down(pool, futures)
        except KeyboardInterrupt:
            stopping = True
    if stopping:
        # set first, so that the running files stop even if an interrupt cuts the log or the wait short
        stop.set()
        log.warning("Stopping the running files at their next chunk. Press Ctrl-C again to stop waiting for them.")
        _shut_down(pool, futures)


def _shut_down(pool: ThreadPoolExecutor, futures: list[Future]) -> None:
    """Start no queued file, and wait for the running files.

    ``pool.shutdown(wait=True)`` would wait as well.
    But in Python 3.12, an interrupt of its ``Thread.join`` marks a running thread as finished,
    so every later join returns at once, including the one at exit.
    """
    pool.shutdown(wait=False, cancel_futures=True)
    # the pool never notifies a future that its shutdown cancelled, so a wait for one would not end
    wait([future for future in futures if not future.cancelled()])
