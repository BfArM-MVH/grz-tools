import os
import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    # For type checking, always use the real signature
    def _sched_getaffinity(pid: int, /) -> set[int]: ...
elif sys.platform == "linux" and hasattr(os, "sched_getaffinity"):
    from os import sched_getaffinity as _sched_getaffinity
else:
    # Fallback for non-Linux platforms
    def _sched_getaffinity(pid: int, /) -> set[int]:
        cpu_count = os.cpu_count() or 1
        return set(range(cpu_count))


def get_cpu_affinity(pid: int = 0) -> set[int]:
    """Return the set of CPUs available to the process (cross-platform safe)."""
    return _sched_getaffinity(pid)


def get_effective_cpu_count(pid: int = 0) -> int:
    """Return the number of CPUs available to the process."""
    return len(get_cpu_affinity(pid))
