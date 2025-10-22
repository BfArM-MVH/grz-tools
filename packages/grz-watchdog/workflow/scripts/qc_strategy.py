import datetime
import math
import random

import yaml
from grz_db.models.submission import SubmissionDb


def get_monthly_submission_stats(db_config_path: str, submitter_id: str) -> tuple[int, int, int]:
    """
    Fetches submission statistics for a given submitter for the current calendar month.
    - Total number of submissions.
    - Number of submissions that have had QC.
    - Number of submissions since the last QC'd submission.
    """
    with open(db_config_path) as f:
        db_url = yaml.safe_load(f)["db"]["database_url"]

    db = SubmissionDb(db_url=db_url, author=None)
    all_submissions = db.list_submissions(limit=None)

    today = datetime.date.today()

    def is_relevant(s):
        return (
            s.submitter_id == submitter_id
            and (d := s.submission_date)
            and d.year == today.year
            and d.month == today.month
        )

    submissions_this_month = sorted(
        list(filter(is_relevant, all_submissions)), key=lambda s: s.submission_date or datetime.date.min
    )

    if not submissions_this_month:
        return 0, 0, 0

    total_this_month = len(submissions_this_month)
    qcd_this_month = sum(1 for sub in submissions_this_month if sub.detailed_qc_passed is not None)

    last_qcd = next(
        filter(lambda item: item[1].detailed_qc_passed is not None, reversed(list(enumerate(submissions_this_month)))),
        None,
    )

    if last_qcd:
        last_qc_index = last_qcd[0]
        submissions_since_last_qc = max(0, total_this_month - 1 - last_qc_index)
    else:
        submissions_since_last_qc = total_this_month

    return total_this_month, qcd_this_month, submissions_since_last_qc


def should_run_qc(db_config_path: str, submitter_id: str, target_percentage: float) -> bool:
    _total_this_month, qcd_this_month, qcd_since_last = get_monthly_submission_stats(db_config_path, submitter_id)

    def _should_run_qc(qcd: int, qcd_since_last: int, target_fraction: float) -> bool:
        if not (0 < target_fraction <= 1):
            raise ValueError("target_fraction must be larger than 0 and at most 1")

        if qcd == 0:
            return True
        block = math.floor(1 / target_fraction)
        n = qcd_since_last + 1
        probability = n / block

        return random.random() < probability

    return _should_run_qc(qcd_this_month, qcd_since_last, target_percentage / 100.0)
