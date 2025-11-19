# packages/grz-watchdog/tests/test_watchdog_interrupts.py

import queue
import subprocess
import tempfile
import threading
import time
from pathlib import Path

import pytest
import yaml

from .conftest import (
    BUCKET_NONCONSENTED,
    INBOX,
    SUBMITTER_ID,
    BaseTest,
    _create_variant_submission,
    run_in_container,
    CONTAINER_RUNTIME,
    GRZ_WATCHDOG_CONTAINER_NAME,
    PIXI_RUN_PREFIX,
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
        print(f"\nRunning to final target: {final_target}, but only until rule decrypt")
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
        print("Intermediate target successfully created.")

        # try to continue from intermediate target
        print(f"\nResuming run to final target: {final_target}")

        result = self._run_watchdog(final_target, config_overrides=config_overrides)
        snakemake_log_output = result.stderr

        print("\nVerifying that initial steps were skipped")
        # check logs to verify snakemake did not re-run the corresponding rules
        assert "localrule download:" not in snakemake_log_output, "Rule 'download' was unexpectedly re-run."
        assert "localrule decrypt:" not in snakemake_log_output, "Rule 'decrypt' was unexpectedly re-run."
        assert "localcheckpoint validate:" in snakemake_log_output, "Rule 'validate' was not run as expected on resume."

        print("OK: Logs confirm that download and decrypt were skipped, and the workflow resumed from validate.")
        print("OK: Logs confirm download and decrypt were skipped.")

        # … then, verify the submission was fully and correctly processed
        self._verify_db_state(submission_id, expected_state="Finished")
        self._verify_archived(submission_id, bucket=BUCKET_NONCONSENTED)

    def test_resumability_after_qc_interrupt(self, test_data_dir: Path, tmp_path: Path, project_root: Path):
        """
        Tests that a workflow interrupted during QC resumes from QC, not from download.
        """
        test_data_dir = Path(test_data_dir) / "wgs_lr"
        submission_dir, submission_id = _create_variant_submission(test_data_dir, "interrupt-qc-test", tmp_path)
        self._submit_data(submission_dir)

        # to have a more reliable interrupt window, try to delay nextflow scripts by inserting `sleep` pre-script execution
        slow_qc_params_host = project_root / "packages/grz-watchdog/tests/config/slow_qc.nextflow_params.yaml"
        slow_qc_params_container = "/tmp/slow_qc.nextflow_params.yaml"
        subprocess.run(
            [
                CONTAINER_RUNTIME,
                "cp",
                str(slow_qc_params_host),
                f"{GRZ_WATCHDOG_CONTAINER_NAME}:{slow_qc_params_container}",
            ],
            check=True,
        )

        config_overrides = {
            "auto-cleanup": "inbox",
            "qc": {
                "selection_strategy": {"enabled": True, "target_percentage": 100.0},
                "run-qc": {"extra": f"-params-file {slow_qc_params_container}"},
            },
        }

        final_target = f"results/{SUBMITTER_ID}/{INBOX}/{submission_id}/processed"
        snakemake_cmd = self._build_snakemake_cmd(final_target, cores=2, config_overrides=config_overrides)

        process_to_interrupt = None
        log_thread = None
        stop_log_thread = threading.Event()

        try:
            # run workflow until we see QC is running, then interrupt
            print(f"\nStarting workflow for {submission_id}, will interrupt during QC")
            process_to_interrupt = self.start_background_process(snakemake_cmd)
            log_queue = queue.Queue()
            log_thread = threading.Thread(
                target=log_consumer, args=(process_to_interrupt.stdout, log_queue, stop_log_thread)
            )
            log_thread.start()

            qc_started = False
            timeout = 180
            start_time = time.time()
            while time.time() - start_time < timeout:
                try:
                    line = log_queue.get(timeout=2)
                    print(line, end="")
                    if "rule qc:" in line:
                        print("\nQC started. Trying forceful interrupt.")
                        self.stop_background_process(process_to_interrupt, force=True)
                        qc_started = True
                        break
                except queue.Empty:
                    if process_to_interrupt.poll() is not None:
                        while not log_queue.empty():
                            print(log_queue.get_nowait(), end="")
                        pytest.fail("Snakemake process exited prematurely before QC started.")

            assert qc_started, f"QC did not start within {timeout} seconds."

            # make sure post-interrupt state is as expected
            print("\nVerifying state after interruption")
            decrypted_dir = f"/workdir/results/{SUBMITTER_ID}/{INBOX}/{submission_id}/decrypted"
            run_in_container("test", "-d", decrypted_dir)
            print("OK: Decrypted directory exists, temp files were preserved.")
            try:
                self._verify_db_state(submission_id, expected_state="QCing")
                print("OK: Database state is 'QCing'.")
            except AssertionError:
                self._verify_db_state(submission_id, expected_state="Error")
                print("OK: Database state is 'Error' after interruption.")

            # after interruption, to unlock/reset snakemake and nextflow
            print("\nUnlocking workflow directory")
            unlock_cmd = self._build_snakemake_cmd(final_target, config_overrides=config_overrides, extra=["--unlock"])
            run_in_container(*unlock_cmd)
            print("OK: Directory unlocked.")

            print("\nCleaning up stale Nextflow work directory")
            nextflow_work_dir = f"/workdir/results/{SUBMITTER_ID}/{INBOX}/{submission_id}/qc/work"
            run_in_container("rm", "-rf", nextflow_work_dir)
            print(f"OK: Removed {nextflow_work_dir}.")

            print("\nResetting DB state to 'reported' to allow resumption")
            run_in_container(
                *PIXI_RUN_PREFIX,
                "grzctl",
                "db",
                "--config-file",
                "/workdir/config/configs/db.yaml",
                "submission",
                "update",
                "--ignore-error-state",
                submission_id,
                "reported",
            )
            print("OK: DB state reset.")

            print(f"\nResuming workflow for {submission_id}")
            resume_config_overrides = config_overrides.copy()
            resume_config_overrides["qc"]["run-qc"].pop("extra", None)

            resume_result = self._run_watchdog(final_target, config_overrides=resume_config_overrides)
            resume_log_output = resume_result.stderr

            # now, verify that we resumed successfully and did not re-do any unnecessary work
            print("\nVerifying logs from resumed run")
            for rule in ["download", "decrypt", "validate", "re_encrypt", "archive", "pruefbericht"]:
                assert f"rule {rule}:" not in resume_log_output
            assert "rule qc:" in resume_log_output
            print("OK: Logs confirm workflow resumed correctly from the QC step.")

            # as usual, verify final outcome as well
            print("\nVerifying final state of submission")
            self._verify_db_state(submission_id, expected_state="Finished")
            self._verify_qc_results_populated(submission_id)
            self._verify_inbox_cleaned(submission_id)
            self._verify_archived(submission_id, bucket=BUCKET_NONCONSENTED)

        finally:
            if process_to_interrupt and process_to_interrupt.poll() is None:
                self.stop_background_process(process_to_interrupt)
            if log_thread:
                stop_log_thread.set()
                log_thread.join(timeout=5)

    def test_resumability_after_pruefbericht_failure(self, test_data_dir: Path, tmp_path: Path):
        """
        Tests that a workflow failing while reporting pruefbericht resumes correctly from that step.
        """
        submission_dir, submission_id = _create_variant_submission(
            test_data_dir / "panel", "failed-pruefbericht-test", tmp_path
        )
        self._submit_data(submission_dir)

        final_target = f"results/{SUBMITTER_ID}/{INBOX}/{submission_id}/processed"

        bad_config_content = yaml.dump(
            {
                "pruefbericht": {
                    "api_base_url": "https://invalid-url.local",
                    "authorization_url": "https://invalid-url.local/token",
                    "client_id": "mock-client-id",
                    "client_secret": "mock-client-password",
                }
            }
        )

        local_bad_config_path = tmp_path / f"{submission_id}_bad_pruefbericht.yaml"
        local_bad_config_path.write_text(bad_config_content)

        container_bad_config_path = f"/tmp/{local_bad_config_path.name}"
        subprocess.run(
            [
                CONTAINER_RUNTIME,
                "cp",
                str(local_bad_config_path),
                f"{GRZ_WATCHDOG_CONTAINER_NAME}:{container_bad_config_path}",
            ],
            check=True,
        )
        config_overrides_fail = {
            "qc": {"selection_strategy": {"enabled": False}},
            "config_paths": {"pruefbericht": container_bad_config_path},
        }

        print(f"\nRunning workflow for {submission_id}, expecting failure at submit_pruefbericht")
        self._run_watchdog_expect_fail(final_target, config_overrides=config_overrides_fail)

        # verify state after the expected failure
        print("\nVerifying state after failure")
        self._verify_db_state(submission_id, expected_state="Error")
        self._verify_archived(submission_id, bucket=BUCKET_NONCONSENTED)
        print("OK: DB state is 'Error' and data is archived, as expected.")

        print("\nUnlocking workflow directory")
        unlock_cmd = self._build_snakemake_cmd(final_target, config_overrides=config_overrides_fail, extra=["--unlock"])
        run_in_container(*unlock_cmd)
        print("OK: Directory unlocked.")

        # reset DB state to "archived" (the one right before "reported")
        print("\nResetting DB state to 'archived' to allow resumption")
        run_in_container(
            *PIXI_RUN_PREFIX,
            "grzctl",
            "db",
            "--config-file",
            "/workdir/config/configs/db.yaml",
            "submission",
            "update",
            "--ignore-error-state",
            submission_id,
            "archived",
        )
        print("OK: DB state reset.")

        # resume the workflow using the default (i.e., correct) pruefbericht config
        print(f"\nResuming workflow for {submission_id}")
        config_overrides_resume = {"qc": {"selection_strategy": {"enabled": False}}}
        resume_result = self._run_watchdog(final_target, config_overrides=config_overrides_resume)
        resume_log_output = resume_result.stderr

        # verify resumption did not unnecessarily re-run previous rules
        print("\nVerifying logs from resumed run")
        for rule in ["download", "decrypt", "validate", "re_encrypt", "archive"]:
            assert f"rule {rule}:" not in resume_log_output, f"Rule '{rule}' was unexpectedly re-run."

        assert "rule submit_pruefbericht:" in resume_log_output, (
            "Rule 'submit_pruefbericht' was not re-run as expected."
        )
        print("OK: Logs confirm workflow resumed correctly from the submit_pruefbericht step.")

        # verify final state
        print("\nVerifying final state of submission")
        self._verify_db_state(submission_id, expected_state="Finished")
        self._verify_inbox_cleaned(submission_id)
        print("OK: Submission successfully processed and cleaned after resumption.")

    def test_resumability_after_omit_from(self, test_data_dir: Path, tmp_path: Path):
        """
        Tests that a workflow run with `--omit-from submit_pruefbericht --notemp` resumes correctly.
        1. Runs the workflow until the 'generate_pruefbericht' step (i.e., omitting from 'submit_pruefbericht' onwards)
        2. Re-runs the workflow to completion.
        3. Verifies that all upstream rules were skipped on the second run.
        """
        submission_dir, submission_id = _create_variant_submission(
            test_data_dir / "panel", "until-pruefbericht-test", tmp_path
        )
        self._submit_data(submission_dir)

        config_overrides = {
            "qc": {"selection_strategy": {"enabled": False}},
        }
        final_target = f"results/{SUBMITTER_ID}/{INBOX}/{submission_id}/processed"
        intermediate_target_file = f"results/{SUBMITTER_ID}/{INBOX}/{submission_id}/pruefbericht.json"

        print(f"\nRunning to final target, but only until rule 'generate_pruefbericht'")
        result = self._run_watchdog(
            final_target,
            config_overrides=config_overrides,
            extra=["--omit-from", "submit_pruefbericht", "--notemp"],
        )
        snakemake_log_output = result.stderr

        # verify which rules ran and which did not
        print("\nVerifying state of partial run")
        for rule in [
            "download",
            "decrypt",
            "re_encrypt",
            "archive",
            "generate_pruefbericht",
        ]:  # "validate" is a checkpoint
            assert f"rule {rule}:" in snakemake_log_output, (
                f"Rule '{rule}' was expected to run but did not. {snakemake_log_output}"
            )
        for rule in ["submit_pruefbericht", "clean", "finalize_success"]:
            assert f"rule {rule}:" not in snakemake_log_output, (
                f"Rule '{rule}' ran but should have been stopped by --until."
            )

        # verify that the intermediate file exists and the DB is in the correct state
        run_in_container("test", "-f", f"/workdir/{intermediate_target_file}")
        self._verify_db_state(submission_id, expected_state="archived")
        print("OK: Intermediate target created and DB state is 'archived'.")

        print(f"\nResuming run to final target: {final_target}")
        resume_result = self._run_watchdog(final_target, config_overrides=config_overrides)
        resume_log_output = resume_result.stderr

        # verify that we did not re-run any steps that were already completed in the previous run
        print("\nVerifying logs from resumed run")
        for rule in ["download", "decrypt", "validate", "re_encrypt", "archive", "generate_pruefbericht"]:
            assert f"rule {rule}:" not in resume_log_output, f"Rule '{rule}' was unexpectedly re-run."

        assert "rule submit_pruefbericht:" in resume_log_output, (
            "Rule 'submit_pruefbericht' was not run on resume as expected."
        )
        assert "rule clean:" in resume_log_output, "Rule 'clean' was not run on resume as expected."
        print("OK: Logs confirm workflow resumed correctly and skipped upstream rules.")

        # verify the final state
        print("\nVerifying final state of submission")
        self._verify_db_state(submission_id, expected_state="Finished")
        self._verify_inbox_cleaned(submission_id)
        print("OK: Submission successfully processed and cleaned after resumption.")
