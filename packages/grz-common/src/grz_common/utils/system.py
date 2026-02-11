import os
import sys

if sys.platform == "linux" and hasattr(os, "sched_getaffinity"):
    from os import sched_getaffinity as _sched_getaffinity
else:

    def _sched_getaffinity(pid=0):
        cpu_count = os.cpu_count() or 1
        return set(range(cpu_count))


def get_cpu_affinity(pid: int = 0) -> set[int]:
    """Return the set of CPUs available to the process (cross-platform safe)."""
    return _sched_getaffinity(pid)
