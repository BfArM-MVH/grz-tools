# packages/grz-watchdog/tests/test_watchdog_interrupts.py

import queue
import subprocess
import threading
import time
from pathlib import Path

import pytest

from .conftest import (
    BUCKET_NONCONSENTED,
    INBOX,
    SUBMITTER_ID,
    BaseTest,
    _create_variant_submission,
    run_in_container,
)


def log_consumer(process_stdout, log_queue, stop_event):
    """Utility to consume log lines from a process in a separate thread."""
    try:
        for line in iter(process_stdout.readline, ""):
            if stop_event.is_set():
                break
            log_queue.put(line)
    finally:
        process_stdout.close()


@pytest.mark.usefixtures("container_test_env")
class TestWorkflowResumption(BaseTest):
    """Tests for workflow resumability after interruptions."""

    def test_resumability_after_decrypt_interrupt(self, test_data_dir: Path, tmp_path: Path):
        """
        Tests that an interrupted workflow resumes correctly.
        1. Runs the workflow up to the decryption step.
        2. Re-runs the workflow to completion.
        3. Verifies that the download and decrypt steps were skipped on the second run.
        """
        test_data_dir = Path(test_data_dir) / "panel"
        submission_dir, submission_id = _create_variant_submission(test_data_dir, "resume-test", tmp_path)
        self._submit_data(submission_dir)

        # to ensure the target directory exists and can be inspected: set auto-cleanup to none
        config_overrides = {
            "auto-cleanup": "none",
            "qc": {"selection_strategy": {"enabled": False}},
        }

        # run to an arbitrary intermediate target by supplying `--until RULE`
        final_target = f"results/{SUBMITTER_ID}/{INBOX}/{submission_id}/processed"
        print(f"\n--- Running to final target: {final_target}, but only until rule decrypt ---")
        result = self._run_watchdog(final_target, config_overrides=config_overrides, extra=["--until", "decrypt"])
        snakemake_log_output = result.stderr

        # check logs to verify snakemake did indeed trigger running some necessary rules
        assert "localrule download:" in snakemake_log_output, "Rule 'download' was not run as expected."
        assert "localrule decrypt:" in snakemake_log_output, "Rule 'decrypt' was not run as expected."
        # but not ones past the target
        assert "localcheckpoint validate:" not in snakemake_log_output, (
            "Rule 'validate' was run, even though it should not have been."
        )

        # verify intermediate dir exists
        intermediate_target = f"results/{SUBMITTER_ID}/{INBOX}/{submission_id}/decrypted"
        run_in_container("test", "-d", f"/workdir/{intermediate_target}")
        print("--- Intermediate target successfully created. ---")

        # try to continue from intermediate target
        print(f"\n--- Resuming run to final target: {final_target} ---")

        result = self._run_watchdog(final_target, config_overrides=config_overrides)
        snakemake_log_output = result.stderr

        print("\n--- Verifying that initial steps were skipped ---")
        # check logs to verify snakemake did not re-run the corresponding rules
        assert "localrule download:" not in snakemake_log_output, "Rule 'download' was unexpectedly re-run."
        assert "localrule decrypt:" not in snakemake_log_output, "Rule 'decrypt' was unexpectedly re-run."
        assert "localcheckpoint validate:" in snakemake_log_output, "Rule 'validate' was not run as expected on resume."

        print("OK: Logs confirm that download and decrypt were skipped, and the workflow resumed from validate.")
        print("OK: Logs confirm download and decrypt were skipped.")

        # … then, verify the submission was fully and correctly processed
        self._verify_db_state(submission_id, expected_state="Finished")
        self._verify_archived(submission_id, bucket=BUCKET_NONCONSENTED)


    def test_resumability_after_qc_failure(self, test_data_dir: Path, tmp_path: Path):
        """
        Tests that a workflow failing during the QC step resumes from QC, not from download.
        1. Runs workflow until decrypted files are available.
        2. Corrupts a FASTQ file to ensure the QC step will fail.
        3. Runs the rest of the workflow and confirms it fails at QC.
        4. Restores the original FASTQ file.
        5. Re-runs the workflow to completion.
        6. Verifies that 'download' and 'decrypt' were skipped and 'qc' was re-run.
        """
        test_data_dir = Path(test_data_dir) / "panel"
        submission_dir, submission_id = _create_variant_submission(test_data_dir, "resume-qc-fail", tmp_path)
        self._submit_data(submission_dir)

        # force QC to run and prevent auto-cleanup of temp files on interruption/failure.
        config_overrides = {
            "auto-cleanup": "none",
            "qc": {"selection_strategy": {"enabled": True, "target_percentage": 100.0}},
        }

        # run until validated
        final_target = f"results/{SUBMITTER_ID}/{INBOX}/{submission_id}/processed"
        self._run_watchdog(final_target, config_overrides=config_overrides, extra=["--until", "validate"])

        decrypted_files_dir = f"/workdir/results/{SUBMITTER_ID}/{INBOX}/{submission_id}/decrypted/files"
        file_to_corrupt = f"{decrypted_files_dir}/S1_R1.fastq.gz"
        backup_file = f"/tmp/{submission_id}_backup.fastq.gz"

        # corrupt one of the decrypted files
        run_in_container("cp", file_to_corrupt, backup_file)
        run_in_container("sh", "-c", f"echo 'garbage' > {file_to_corrupt}")

        # resume run, should fail
        self._run_watchdog_expect_fail(final_target, config_overrides=config_overrides)
        self._verify_db_state(submission_id, expected_state="Error")

        # reinstate backuped file, re-run (update db state first to be able to continue)
        run_in_container("cp", backup_file, file_to_corrupt)
        run_in_container(
            *self._build_snakemake_cmd(final_target)[1:-1],
            "grzctl", "db", "--config-file", "/workdir/config/configs/db.yaml",
            "submission", "update", submission_id, "validated"
        )
        resume_result = self._run_watchdog(final_target, config_overrides=config_overrides)
        resume_log_output = resume_result.stderr

        assert "localrule download:" not in resume_log_output, "Rule 'download' was unexpectedly re-run."
        assert "localrule decrypt:" not in resume_log_output, "Rule 'decrypt' was unexpectedly re-run."
        assert "rule qc:" in resume_log_output, "Rule 'qc' was not re-run as expected."

        self._verify_db_state(submission_id, expected_state="Finished")
        self._verify_qc_results_populated(submission_id)
        self._verify_inbox_cleaned(submission_id)
        self._verify_archived(submission_id, bucket=BUCKET_NONCONSENTED)