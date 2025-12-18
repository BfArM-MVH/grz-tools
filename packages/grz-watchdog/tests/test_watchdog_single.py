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


@pytest.mark.usefixtures("container_test_env")
class TestProcessSingle(BaseTest):
    """Testing setup for the process-single branch of grz-watchdog"""

    @pytest.mark.parametrize("submission_type, submission_id", TEST_CASES)
    def test_single_valid_submission(
        self, test_data_dir: Path, tmp_path: Path, submission_type: str, submission_id: str
    ):
        """
        Test the successful end-to-end processing of a single valid submission.
        """
        test_data_dir = Path(test_data_dir) / submission_type
        valid_dir, submission_id = _create_variant_submission(test_data_dir, "initial", tmp_path)
        self._submit_data(valid_dir)

        target_file = f"results/{SUBMITTER_ID}/{INBOX}/{submission_id}/processed"
        config_overrides = {"qc": {"selection_strategy": {"enabled": False}}}
        self._run_watchdog(target_file, config_overrides=config_overrides)

        self._verify_db_state(submission_id, expected_state="Finished")
        self._verify_inbox_cleaned(submission_id)
        self._verify_archived(submission_id, bucket=BUCKET_NONCONSENTED)

    def test_single_invalid_submission(self, test_data_dir: Path, tmp_path: Path):
        """
        Test failure path for submission with corrupted (encrypted) fastq file.
        """
        test_data_dir = Path(test_data_dir) / "panel"
        variant_corrupted_dir, submission_id = _create_variant_submission(test_data_dir, "corrupted", tmp_path)
        self._submit_data(variant_corrupted_dir)

        metadata_path = variant_corrupted_dir / "metadata" / "metadata.json"
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

        run_in_container(*PIXI_RUN_PREFIX, "mc", "cp", s3_path_to_corrupt, container_temp_path)
        corruption_command = f"printf '\\1' | dd of={container_temp_path} bs=1 count=1 conv=notrunc"
        run_in_container("sh", "-c", corruption_command)
        run_in_container(*PIXI_RUN_PREFIX, "mc", "cp", container_temp_path, s3_path_to_corrupt)

        target_file = f"results/{SUBMITTER_ID}/{INBOX}/{submission_id}/processed"
        config_overrides = {"qc": {"selection_strategy": {"enabled": False}}}
        self._run_watchdog_expect_fail(target_file, config_overrides=config_overrides)

        self._verify_db_state(submission_id, expected_state="Error")

        try:
            result = run_in_container(*PIXI_RUN_PREFIX, "mc", "ls", f"adm/{BUCKET_INBOX}/{submission_id}")
            assert any("files/" in l for l in result.stdout.strip().splitlines())
        except subprocess.CalledProcessError as e:
            pytest.fail(f"Failed to check inbox state for {submission_id}. Error: {e.stderr}")

    def test_single_valid_submission_with_qc(self, test_data_dir: Path, tmp_path: Path):
        """
        Test the successful end-to-end processing of a single valid submission including the QC pipeline.
        """
        submission_type = "panel"
        test_data_dir = Path(test_data_dir) / submission_type
        variant_qc_dir, submission_id = _create_variant_submission(test_data_dir, "qc", tmp_path)
        self._submit_data(variant_qc_dir)

        target_file = f"results/{SUBMITTER_ID}/{INBOX}/{submission_id}/processed"
        config_overrides = {"qc": {"selection_strategy": {"enabled": True, "target_percentage": 100.0}}}
        self._run_watchdog(target_file, config_overrides=config_overrides)

        self._verify_db_state(submission_id, expected_state="Finished")
        self._verify_qc_results_populated(submission_id)
        self._verify_inbox_cleaned(submission_id)
        self._verify_archived(submission_id, bucket=BUCKET_NONCONSENTED)

    def test_valid_submission_after_failed_submission(self, test_data_dir: Path, tmp_path: Path):
        """
        Tests that the workflow can process a new, valid submission even after a
        previous run failed, leaving potentially stale intermediate files.
        """
        config_overrides = {"qc": {"selection_strategy": {"enabled": False}}}
        base_submission_dir = test_data_dir / "panel"
        fail_dir, fail_id = _create_variant_submission(base_submission_dir, "fail_run", tmp_path)
        self._submit_data(fail_dir)

        metadata_path = fail_dir / "metadata" / "metadata.json"
        with open(metadata_path) as f:
            metadata = GrzSubmissionMetadata(**json.load(f))

        file_to_corrupt_meta = next(
            file_meta
            for donor in metadata.donors
            for lab_datum in donor.lab_data
            if lab_datum.sequence_data
            for file_meta in lab_datum.sequence_data.files
            if file_meta.file_type == FileType.fastq
        )

        s3_path_to_corrupt = f"adm/{BUCKET_INBOX}/{fail_id}/files/{file_to_corrupt_meta.encrypted_file_path()}"
        container_temp_path = f"/tmp/{Path(file_to_corrupt_meta.file_path).name}.corrupt"
        run_in_container(*PIXI_RUN_PREFIX, "mc", "cp", s3_path_to_corrupt, container_temp_path)
        run_in_container("sh", "-c", f"printf '\\0' | dd of={container_temp_path} bs=1 count=1 conv=notrunc")
        run_in_container(*PIXI_RUN_PREFIX, "mc", "cp", container_temp_path, s3_path_to_corrupt)

        fail_target = f"results/{SUBMITTER_ID}/{INBOX}/{fail_id}/processed"
        self._run_watchdog_expect_fail(fail_target, config_overrides=config_overrides)
        self._verify_db_state(fail_id, expected_state="Error")

        pass_dir, pass_id = _create_variant_submission(base_submission_dir, "pass_run", tmp_path)
        self._submit_data(pass_dir)

        pass_target = f"results/{SUBMITTER_ID}/{INBOX}/{pass_id}/processed"
        self._run_watchdog(pass_target, config_overrides=config_overrides)

        self._verify_db_state(pass_id, expected_state="Finished")
        self._verify_inbox_cleaned(pass_id)
        self._verify_archived(pass_id, bucket=BUCKET_NONCONSENTED)

    def test_single_submission_with_existing_db(self, test_data_dir: Path, tmp_path: Path):
        """
        Tests that processing a new submission does not affect submissions already
        present and finished in the database.
        """
        base_submission_dir = test_data_dir / "panel"
        config_overrides = {"qc": {"selection_strategy": {"enabled": False}}}

        pre_existing_dir, pre_existing_id = _create_variant_submission(base_submission_dir, "pre_existing", tmp_path)
        self._submit_data(pre_existing_dir)

        pre_existing_target = f"results/{SUBMITTER_ID}/{INBOX}/{pre_existing_id}/processed"
        self._run_watchdog(pre_existing_target, config_overrides=config_overrides)

        self._verify_db_state(pre_existing_id, expected_state="Finished")
        self._verify_archived(pre_existing_id, bucket=BUCKET_NONCONSENTED)
        self._verify_inbox_cleaned(pre_existing_id)

        new_submission_dir, new_submission_id = _create_variant_submission(
            base_submission_dir, "new_submission", tmp_path
        )
        self._submit_data(new_submission_dir)

        new_target = f"results/{SUBMITTER_ID}/{INBOX}/{new_submission_id}/processed"
        self._run_watchdog(new_target, config_overrides=config_overrides)

        self._verify_db_state(new_submission_id, expected_state="Finished")
        self._verify_archived(new_submission_id, bucket=BUCKET_NONCONSENTED)
        self._verify_inbox_cleaned(new_submission_id)

        self._verify_db_state(pre_existing_id, expected_state="Finished")
