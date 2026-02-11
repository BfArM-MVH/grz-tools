import threading
from collections import defaultdict
from typing import Any


class SubmissionContext:
    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._stats: dict[str, dict[str, Any]] = defaultdict(dict)
        self._errors: list[str] = []

    def record_stats(self, file_path: str, stats: dict[str, Any]) -> None:
        with self._lock:
            self._stats[file_path].update(stats)

    def get_stats(self, file_path: str) -> dict[str, Any]:
        with self._lock:
            return self._stats.get(file_path, {}).copy()

    def add_error(self, error: str) -> None:
        with self._lock:
            self._errors.append(error)

    @property
    def has_errors(self) -> bool:
        with self._lock:
            return len(self._errors) > 0


class ConsistencyValidator:
    def __init__(self, context: SubmissionContext):
        self.context = context

    def check_pair(self, path_a: str, path_b: str) -> bool:
        """
        Fail-Fast check: If both files are done, compare them.
        Returns False if mismatch.
        """
        stats_a = self.context.get_stats(path_a)
        stats_b = self.context.get_stats(path_b)

        # If one is missing, it's not ready to check (or failed earlier). Pass for now.
        if not stats_a or not stats_b:
            return True

        reads_a = stats_a.get("read_count")
        reads_b = stats_b.get("read_count")

        # Only compare if both successfully counted reads
        if reads_a is not None and reads_b is not None and reads_a != reads_b:
            msg = f"Read Count Mismatch: {path_a} ({reads_a}) != {path_b} ({reads_b})"
            self.context.add_error(msg)
            return False
        return True
