import subprocess
from pathlib import Path

import pytest

from .conftest import (
    BUCKET_CONSENTED,
    BUCKET_INBOX,
    BUCKET_NONCONSENTED,
    CONTAINER_COMPOSE_CMD,
    CONTAINER_RUNTIME,
    DOCKER_COMPOSE_FILE,
    GRZ_SUBMITTER_CONTAINER_NAME,
    GRZ_SUBMITTER_SERVICE_NAME,
    GRZ_WATCHDOG_CONTAINER_NAME,
    GRZ_WATCHDOG_SERVICE_NAME,
    INBOX,
    MINIO_SERVICE_NAME,
    PIXI_RUN_PREFIX,
    SNAKEMAKE_BASE_CMD,
    SUBMITTER_ID,
    run_in_container,
    BaseTest,
)

TEST_CASES = [
    ("panel", "123456789_2024-11-08_d0f805c5"),
    # ("wgs", "123456789_2024-10-28_e1bab61b"),
    ("wgs_lr", "123456789_2024-10-28_e1bab61b"),
]


@pytest.mark.usefixtures("setup_class_environment")
class TestProcessSingle(BaseTest):
    """Testing setup for the process-single branch of grz-watchdog"""

    def _submit_data(self, test_data_dir: Path, submission_type: str):
        """Helper to prepare and submit a specific test dataset."""
        print(f"\n--- Submitting Data for test case: {submission_type} ---")
        local_data_path = test_data_dir / submission_type
        container_data_path = f"/tmp/{submission_type}"

        subprocess.run(
            [CONTAINER_RUNTIME, "cp", str(local_data_path), f"{GRZ_SUBMITTER_CONTAINER_NAME}:{container_data_path}"],
            check=True,
        )
        run_in_container(
            *PIXI_RUN_PREFIX,
            "grz-cli",
            "submit",
            "--submission-dir",
            container_data_path,
            "--config-file",
            "/workdir/config/grz-cli.config.yaml",
            service=GRZ_SUBMITTER_SERVICE_NAME,
        )

    @pytest.mark.parametrize("submission_type, submission_id", TEST_CASES)
    def test_single_valid_submission(self, test_data_dir: Path, submission_type: str, submission_id: str):
        """
        Tests the successful end-to-end processing of a single valid submission.
        """
        self._submit_data(test_data_dir, submission_type)

        target_file = f"results/{SUBMITTER_ID}/{INBOX}/{submission_id}/processed/without_qc"
        self._run_snakemake(target_file)

        self._verify_db_state(submission_id, expected_state="Finished")
        self._verify_inbox_cleaned(submission_id)
        self._verify_archived(submission_id, bucket=BUCKET_CONSENTED)
