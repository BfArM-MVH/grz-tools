import json
import subprocess
from pathlib import Path

import pytest
from grz_pydantic_models.submission.metadata import FileType, GrzSubmissionMetadata

from .conftest import (
    BUCKET_INBOX,
    BUCKET_NONCONSENTED,
    INBOX,
    PIXI_RUN_PREFIX,
    SUBMITTER_ID,
    BaseTest,
    _create_variant_submission,
    run_in_container,
)

TEST_CASES = [
    ("panel", "123456789_2024-11-08_d0f805c5"),
    # ("wgs", "123456789_2024-10-28_e1bab61b"),
    ("wgs_lr", "123456789_2024-10-28_e1bab61b"),
]


@pytest.mark.usefixtures("setup_class_environment")
class TestProcessSingle(BaseTest):
    """Testing setup for the process-single branch of grz-watchdog"""

    @pytest.mark.parametrize("submission_type, submission_id", TEST_CASES)
    def test_single_valid_submission(self, test_data_dir: Path, submission_type: str, submission_id: str):
        """
        Test the successful end-to-end processing of a single valid submission.
        """
        self._submit_data(test_data_dir / submission_type)

        target_file = f"results/{SUBMITTER_ID}/{INBOX}/{submission_id}/processed/without_qc"
        self._run_watchdog(target_file)

        self._verify_db_state(submission_id, expected_state="Finished")
        self._verify_inbox_cleaned(submission_id)
        self._verify_archived(submission_id, bucket=BUCKET_NONCONSENTED)

    def test_single_invalid_submission(self, test_data_dir: Path, tmp_path: Path):
        """
        Test failure path for submission with corrupted (encrypted) fastq file.
        """
        test_data_dir = Path(test_data_dir) / "panel"
        variant3_dir, sub_id3 = _create_variant_submission(test_data_dir, "corrupted", tmp_path)
        submission_id = sub_id3

        self._submit_data(variant3_dir)

        metadata_path = variant3_dir / "metadata" / "metadata.json"
        with open(metadata_path) as f:
            metadata = GrzSubmissionMetadata(**json.load(f))

        file_to_corrupt_meta = None
        for donor in metadata.donors:
            for lab_datum in donor.lab_data:
                if lab_datum.sequence_data:
                    for file_meta in lab_datum.sequence_data.files:
                        if file_meta.file_type == FileType.fastq:
                            file_to_corrupt_meta = file_meta
                            break
                    if file_to_corrupt_meta:
                        break
            if file_to_corrupt_meta:
                break

        assert file_to_corrupt_meta is not None, "Could not find a FASTQ file to corrupt in the metadata."

        encrypted_file_key = f"files/{file_to_corrupt_meta.encrypted_file_path()}"

        s3_path_to_corrupt = f"adm/{BUCKET_INBOX}/{submission_id}/{encrypted_file_key}"
        container_temp_path = f"/tmp/{Path(encrypted_file_key).name}"

        run_in_container(
            *PIXI_RUN_PREFIX,
            "mc",
            "alias",
            "set",
            "adm",
            "http://minio:9000",
            "minioadmin",
            "minioadmin",
        )
        run_in_container(*PIXI_RUN_PREFIX, "mc", "cp", s3_path_to_corrupt, container_temp_path)
        corruption_command = f"printf '\\1' | dd of={container_temp_path} bs=1 count=1 conv=notrunc"
        run_in_container("sh", "-c", corruption_command)
        run_in_container(*PIXI_RUN_PREFIX, "mc", "cp", container_temp_path, s3_path_to_corrupt)

        target_file = f"results/{SUBMITTER_ID}/{INBOX}/{submission_id}/processed/without_qc"
        self._run_watchdog_expect_fail(target_file)

        self._verify_db_state(submission_id, expected_state="Error")

        try:
            result = run_in_container(*PIXI_RUN_PREFIX, "mc", "ls", f"adm/{BUCKET_INBOX}/{submission_id}")
            assert any("files/" in l for l in result.stdout.strip().splitlines())
        except subprocess.CalledProcessError as e:
            pytest.fail(f"Failed to check inbox state for {submission_id}. Error: {e.stderr}")
