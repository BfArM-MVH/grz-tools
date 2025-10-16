from pathlib import Path

import pytest

from .conftest import (
    BUCKET_NONCONSENTED,
    BaseTest,
    _create_variant_submission,
)


@pytest.mark.usefixtures("setup_class_environment")
class TestProcessMulti(BaseTest):
    """Testing setup for batch, daemon, and failure modes of grz-watchdog"""

    def test_batch_valid_submissions(self, test_data_dir: Path, tmp_path: Path):
        """
        Test successful processing of multiple valid submissions in batch mode.
        """
        test_data_dir, tmp_path = Path(test_data_dir), Path(tmp_path)
        base_submission_dir = test_data_dir / "panel"

        # create and submit two submissions
        variant1_dir, sub_id1 = _create_variant_submission(base_submission_dir, "v1", tmp_path)
        variant2_dir, sub_id2 = _create_variant_submission(base_submission_dir, "v2", tmp_path)
        self._submit_data(variant1_dir)
        self._submit_data(variant2_dir)

        # run watchdog with target "pending", which should process all available submissions
        self._run_watchdog("pending")

        # verify both submissions have been processed
        for sub_id in [sub_id1, sub_id2]:
            self._verify_db_state(sub_id, expected_state="Finished")
            self._verify_inbox_cleaned(sub_id)
            self._verify_archived(sub_id, bucket=BUCKET_NONCONSENTED)
