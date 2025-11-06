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


    def test_resumability_after_qc_interrupt(self, test_data_dir: Path, tmp_path: Path):
        """
        Tests that a workflow interrupted during the QC step resumes from QC, not from download.
        1. Runs a workflow with QC enabled in the background.
        2. Interrupts the process once the 'qc' rule starts.
        3. Verifies that intermediate decrypted files are preserved.
        4. Re-runs the workflow to completion.
        5. Verifies that the 'download' and 'decrypt' steps were skipped on the second run,
           and that 'qc' was attempted again.
        """
        test_data_dir = Path(test_data_dir) / "panel"
        submission_dir, submission_id = _create_variant_submission(test_data_dir, "resume-qc-test", tmp_path)
        self._submit_data(submission_dir)

        # Force QC to run and prevent auto-cleanup of temp files on interruption.
        config_overrides = {
            "auto-cleanup": "none",
            "qc": {"selection_strategy": {"enabled": True, "target_percentage": 100.0}},
        }

        final_target = f"results/{SUBMITTER_ID}/{INBOX}/{submission_id}/processed"
        snakemake_cmd = self._build_snakemake_cmd(final_target, cores=2, config_overrides=config_overrides)

        process_to_interrupt = None
        log_thread = None
        stop_log_thread = threading.Event()

        try:
            # --- Part 1: Run and Interrupt ---
            print(f"\n--- Running workflow for {submission_id} and interrupting during QC ---")
            process_to_interrupt = self.start_background_process(snakemake_cmd)
            log_queue = queue.Queue()
            log_thread = threading.Thread(
                target=log_consumer, args=(process_to_interrupt.stdout, log_queue, stop_log_thread)
            )
            log_thread.start()

            qc_started = False
            start_time = time.time()
            timeout = 180  # Timeout for waiting for QC to start
            while time.time() - start_time < timeout:
                try:
                    line = log_queue.get(timeout=2)
                    print(line, end="")
                    if "rule qc:" in line:
                        print("\n--- QC rule started. Interrupting workflow now. ---")
                        self.stop_background_process(process_to_interrupt)
                        qc_started = True
                        break
                except queue.Empty:
                    if process_to_interrupt.poll() is not None:
                        pytest.fail("Snakemake process exited prematurely before QC started.")

            assert qc_started, f"QC rule did not start within {timeout} seconds."

            # --- Verification after Part 1 ---
            print("\n--- Verifying state after interruption ---")
            decrypted_dir = f"/workdir/results/{SUBMITTER_ID}/{INBOX}/{submission_id}/decrypted"
            run_in_container("test", "-d", decrypted_dir)
            print("OK: Decrypted directory exists, temp files were preserved.")

            qc_output_dir = f"/workdir/results/{SUBMITTER_ID}/{INBOX}/{submission_id}/qc/out"
            with pytest.raises(subprocess.CalledProcessError):
                run_in_container("test", "-d", qc_output_dir)
            print("OK: QC output directory does not exist, confirming interruption.")

            # --- Part 2: Resume Workflow ---
            print(f"\n--- Resuming workflow for {submission_id} ---")
            resume_result = self._run_watchdog(final_target, config_overrides=config_overrides)
            resume_log_output = resume_result.stderr

            # --- Verification after Part 2 ---
            print("\n--- Verifying logs from resumed run ---")
            assert "localrule download:" not in resume_log_output, "Rule 'download' was unexpectedly re-run."
            assert "localrule decrypt:" not in resume_log_output, "Rule 'decrypt' was unexpectedly re-run."
            assert "rule qc:" in resume_log_output, "Rule 'qc' was not re-run as expected on resume."
            print("OK: Logs confirm workflow resumed from the QC step.")

            # --- Final State Verification ---
            print("\n--- Verifying final state of submission ---")
            self._verify_db_state(submission_id, expected_state="Finished")
            self._verify_qc_results_populated(submission_id)
            self._verify_inbox_cleaned(submission_id)
            self._verify_archived(submission_id, bucket=BUCKET_NONCONSENTED)
            print("\n--- Test successful! ---")

        finally:
            if process_to_interrupt:
                self.stop_background_process(process_to_interrupt)
            if log_thread:
                stop_log_thread.set()
                log_thread.join(timeout=5)
