import queue
import sys
import threading
from collections.abc import Callable, Iterable
from concurrent.futures import ThreadPoolExecutor
from typing import Any

from grz_common.constants import TQDM_DEFAULTS
from tqdm.auto import tqdm


def _run_parallel_with_progress[T](
    items: Iterable[T],
    get_size_fn: Callable[[T], int],
    worker_fn: Callable[[T, tqdm, threading.Lock, Any], None],
    threads: int,
    global_desc: str = "GLOBAL  ",
) -> None:
    """Generic runner handling ThreadPool, a global progress bar, and local progress bars."""
    items_list = list(items)
    if not items_list:
        return

    total_bytes = sum(get_size_fn(item) for item in items_list)
    num_workers = min(threads, len(items_list))

    defaults = TQDM_DEFAULTS.copy()
    defaults["leave"] = True

    global_lock = threading.Lock()
    worker_pbars = []
    for i in range(1, num_workers + 1):
        worker_pbars.append(tqdm(total=0, desc="WAITING ", position=i, **defaults))  # type: ignore[call-overload]

    pbar_queue: queue.Queue[tqdm] = queue.Queue()
    for pbar in worker_pbars:
        pbar_queue.put(pbar)

    with (
        tqdm(total=total_bytes, desc=global_desc, position=0, file=sys.stderr, **TQDM_DEFAULTS) as pbar_global,  # type: ignore[call-overload]
        ThreadPoolExecutor(max_workers=threads) as pool,
    ):

        def _thread_wrapper(item: T):
            pbar_local = pbar_queue.get()
            try:
                worker_fn(item, pbar_local, global_lock, pbar_global)
            finally:
                pbar_local.set_postfix({})
                pbar_local.reset(total=0)
                pbar_local.set_description("WAITING ")
                pbar_local.refresh()
                pbar_queue.put(pbar_local)

        futures = [pool.submit(_thread_wrapper, item) for item in items_list]

        for future in futures:
            future.result()

    for pbar in worker_pbars:
        pbar.close()
