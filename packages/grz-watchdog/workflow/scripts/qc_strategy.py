import datetime
import math
import random

import yaml
from grz_db.models.submission import SubmissionDb


def should_run_qc(db_config_path: str, submitter_id: str, target_percentage: float, salt: str) -> bool:
    """
    Determines if a validated submission should undergo QC.
    """
    with open(db_config_path) as f:
        db_url = yaml.safe_load(f)["db"]["database_url"]
    db = SubmissionDb(db_url=db_url, author=None)

    today = datetime.date.today()
    num_validated, num_qcing, num_since_last_qcing = db.get_monthly_qc_stats(
        submitter_id=submitter_id, year=today.year, month=today.month
    )

    # always qc the first (validated) submission each month
    if num_validated == 0:
        return True

    # for example, for a target_percentage of 2%, pick 1 submission at random from 50 submissions
    target_fraction = target_percentage / 100.0
    block = math.floor(1 / target_fraction)
    block_index = num_qcing
    current_index_in_block = num_since_last_qcing

    # ensure the seed is consistent across submitter, date, block_index and (secret) salt
    seed = f"{submitter_id}-{today.year}-{today.month}-{block_index}-{salt}"
    rng = random.Random(seed)  # noqa: S311
    target_index_in_block = rng.randint(0, block - 1)

    if current_index_in_block == target_index_in_block:
        return True

    # if we somehow have exceeded the block size, return True, otherwise False
    return current_index_in_block >= block - 1
