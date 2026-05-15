import queue
import threading
from collections.abc import Callable, Iterable
from concurrent.futures import ThreadPoolExecutor
from typing import Any

from grz_common.constants import TQDM_DEFAULTS
from tqdm.auto import tqdm


def _run_parallel_with_progress[T](
    items: Iterable[T],
    get_size_fn: Callable[[T], int],
    worker_fn: Callable[[T, int, threading.Lock, tqdm | Any], None],
    threads: int,
    global_desc: str = "Global Progress",
) -> None:
    """Generic runner handling a ThreadPool, a global progress bar, and a thread-safe queue for fixed UI positions."""
    items_list = list(items)
    if not items_list:
        return

    total_bytes = sum(get_size_fn(item) for item in items_list)
    num_workers = min(threads, len(items_list))

    position_queue: queue.Queue[int] = queue.Queue()
    for i in range(1, num_workers + 1):
        position_queue.put(i)

    global_lock = threading.Lock()

    with (
        tqdm(total=total_bytes, desc=global_desc, position=0, **TQDM_DEFAULTS) as pbar_global,  # type: ignore[call-overload]
        ThreadPoolExecutor(max_workers=threads) as pool,
    ):

        def _thread_wrapper(item: T):
            pos = position_queue.get()
            try:
                worker_fn(item, pos, global_lock, pbar_global)
            finally:
                position_queue.put(pos)

        futures = [pool.submit(_thread_wrapper, item) for item in items_list]

        for future in futures:
            future.result()
