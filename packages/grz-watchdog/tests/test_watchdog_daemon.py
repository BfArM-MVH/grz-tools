import queue
import threading
import time
from pathlib import Path

import pytest

from .conftest import BUCKET_NONCONSENTED, BaseTest, _create_variant_submission


def log_consumer(process_stdout, log_queue, stop_event):
    try:
        for line in iter(process_stdout.readline, ""):
            if stop_event.is_set():
                break
            log_queue.put(line)
    finally:
        process_stdout.close()

@pytest.mark.usefixtures("container_test_env")
class TestDaemonMode(BaseTest):
    """Testing the continuous monitoring (daemon) mode of grz-watchdog."""

    def test_daemon_picks_up_new_submission(self, test_data_dir: Path, tmp_path: Path):
        """
        Test that the daemon, when running, detects and processes a newly arrived submission.
        """
        config_overrides = {"qc": {"selection_strategy": {"enabled": False}}}
        stop_log_thread = threading.Event()
        log_thread = None
        process = None

        try:
            process = self.start_background_process(
                self._build_snakemake_cmd("daemon", cores=2, config_overrides=config_overrides)
            )
            log_queue = queue.Queue()
            log_thread = threading.Thread(
                target=log_consumer, args=(process.stdout, log_queue, stop_log_thread)
            )
            log_thread.start()

            ready = False
            start_time = time.time()
            timeout = 60

            while time.time() - start_time < timeout:
                try:
                    line = log_queue.get(timeout=1)
                    print(line, end="")
                    if "Starting monitoring thread" in line:
                        ready = True
                        break
                except queue.Empty:
                    if process.poll() is not None:
                        pytest.fail("Snakemake daemon process exited prematurely.")
                    continue

            assert ready, f"Daemon did not become ready within {timeout} seconds."

            base_submission_dir = test_data_dir / "panel"
            submission_dir, submission_id = _create_variant_submission(base_submission_dir, "daemon_test", tmp_path)
            self._submit_data(submission_dir)

            success = False
            start_time = time.time()
            processing_timeout = 180

            while time.time() - start_time < processing_timeout:
                try:
                    while True:
                        line = log_queue.get_nowait()
                        print(line, end="")
                except queue.Empty:
                    pass

                try:
                    self._verify_db_state(submission_id, expected_state="Finished")
                    print(f"--- Submission {submission_id} successfully processed! ---")
                    success = True
                    break
                except (AssertionError, Exception):
                    time.sleep(10)

            assert success, f"Submission {submission_id} was not processed by the daemon within {processing_timeout} seconds."

            self._verify_inbox_cleaned(submission_id)
            self._verify_archived(submission_id, bucket=BUCKET_NONCONSENTED)

        finally:
            if process:
                self.stop_background_process(process)
            if log_thread:
                stop_log_thread.set()
                log_thread.join(timeout=5)