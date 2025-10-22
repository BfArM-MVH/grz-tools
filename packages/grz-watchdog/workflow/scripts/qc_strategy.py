import calendar
import datetime
import math
import random

import yaml
from grz_db.models.submission import SubmissionDb


def get_monthly_submission_stats(db_config_path: str, submitter_id: str) -> tuple[int, int]:
    """
    Fetches submission statistics for a given submitter for the current calendar month.

    Returns:
        A tuple containing:
        - Total number of submissions this month.
        - Number of qc'd submissions this month.
    """
    with open(db_config_path) as f:
        db_url = yaml.safe_load(f)["db"]["database_url"]

    db = SubmissionDb(db_url=db_url, author=None)
    all_submissions = db.list_submissions(limit=None)

    today = datetime.date.today()
    total_this_month = 0
    qc_this_month = 0

    for sub in all_submissions:
        if sub.submitter_id != submitter_id:
            continue

        if (d := sub.submission_date) and d.year == today.year and d.month == today.month:
            total_this_month += 1
            if sub.detailed_qc_passed is not None:
                qc_this_month += 1

    return total_this_month, qc_this_month


def should_run_qc(db_config_path: str, submitter_id: str, target_percentage: float) -> bool:
    total_this_month, qc_this_month = get_monthly_submission_stats(db_config_path, submitter_id)

    def _should_run_qc(total: int, qcd: int, target_fraction: float) -> bool:
        # first submission of the month for an LE always gets qc'd.
        if total == 0:
            return True

        # if no submissions have been qc'd this month, pick the next one.
        if qcd == 0:
            return True

        # increment total by 1, as at this point a submission is currently being added
        total = total + 1
        required_qc_count = max(1, math.ceil(total * target_fraction))

        # if we are already at or above the required count, just "flip a coin" with the base probability.
        if qc_this_month >= required_qc_count:
            return random.random() < target_fraction  # noqa: S311

        # if we are behind, calculate a reasonable probability needed to catch up.
        additional_qcs_needed = required_qc_count - qc_this_month

        # estimate remaining submissions based on current rate and time left in the month.
        today = datetime.date.today()
        _day1, days_in_month = calendar.monthrange(today.year, today.month)
        day_of_month = today.day

        if day_of_month > 0:
            rate = total / day_of_month
            days_remaining = days_in_month - day_of_month
            estimated_remaining_submissions = rate * days_remaining
        else:
            estimated_remaining_submissions = 1

        pseudo_prob = additional_qcs_needed / max(1, estimated_remaining_submissions)  # can be larger than 1
        return random.random() < pseudo_prob  # noqa: S311

    return _should_run_qc(total_this_month, qc_this_month, target_percentage / 100.0)
