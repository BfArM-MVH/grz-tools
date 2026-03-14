import logging
import threading
from collections import defaultdict
from typing import Any, Self

from grz_common.workers.submission import SubmissionMetadata

log = logging.getLogger(__name__)


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


class ReadPairConsistencyValidator:
    def __init__(self, context: SubmissionContext, partner_map: dict[str, str]) -> None:
        self.context = context
        self.partner_map = partner_map

    @staticmethod
    def get_partner_map(submission_metadata: SubmissionMetadata) -> dict[str, str]:
        partner_map = {}
        for _donor, _lab_datum, pairs, _thresholds in submission_metadata.iter_paired_end_fastqs():
            for fq1, fq2 in pairs:
                partner_map[fq1.file_path] = fq2.file_path
                partner_map[fq2.file_path] = fq1.file_path
        return partner_map

    @classmethod
    def from_submission_metadata(cls, context: SubmissionContext, submission_metadata: SubmissionMetadata) -> Self:
        return cls(context, cls.get_partner_map(submission_metadata))

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

    def check(self, path: str) -> bool:
        partner = self.partner_map.get(path)
        if partner is None:
            log.debug(f"Partner for '{path}' not found")
            return True
        return self.check_pair(partner, path)
