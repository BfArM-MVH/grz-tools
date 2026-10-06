"""Stop the files that a pool processes in parallel, one level further on each interrupt."""

import logging
import threading
from collections.abc import Iterable
from concurrent.futures import Future, ThreadPoolExecutor, as_completed

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

    :param pool: The pool that processes the files.
    :param futures: The futures of the files in ``pool``.
    :param stop: The event that the files check before each chunk.
    :raises KeyboardInterrupt: The first interrupt, or the third one.
    """
    try:
        for future in as_completed(futures):
            future.result()
    except KeyboardInterrupt as interrupt:
        stopping = isinstance(interrupt, TerminateInterrupt)
        if not stopping:
            log.warning(
                "Interrupted: no queued file starts, and the running files finish. Press Ctrl-C again to stop them."
            )
            try:
                pool.shutdown(wait=True, cancel_futures=True)
            except KeyboardInterrupt:
                stopping = True
        if stopping:
            log.warning("Stopping the running files at their next chunk. Press Ctrl-C again to stop waiting for them.")
            stop.set()
            pool.shutdown(wait=True, cancel_futures=True)
        raise
    except Exception:
        pool.shutdown(cancel_futures=True)
        raise
