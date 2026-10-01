"""Tests for submission_id-aware validation skip logic"""

import time
from pathlib import Path

import pytest
from grz_common.progress.progress_logging import FileProgressLogger
from grz_common.progress.states import ValidationState
from grz_common.workers.submission import Submission


@pytest.fixture
def submission(submission_metadata_dir) -> Submission:
    return Submission(
        metadata_dir=submission_metadata_dir,
        files_dir=Path("tests/mock_files/submissions/valid_submission/files"),
    )


@pytest.fixture
def temp_checksum_log(tmp_path) -> Path:
    return tmp_path / "progress_checksum.cjson"


@pytest.fixture
def temp_seq_data_log(tmp_path) -> Path:
    return tmp_path / "progress_seq_data.cjson"


def _seed_validation_logs(submission, checksum_log, seq_data_log, submission_id, validation_passed=True):
    """Pre-seed both progress logs for all files in the submission."""
    checksum_logger = FileProgressLogger[ValidationState](log_file_path=checksum_log)
    seq_data_logger = FileProgressLogger[ValidationState](log_file_path=seq_data_log)
    for file_path, file_metadata in submission.files.items():
        state = ValidationState(validation_passed=validation_passed, submission_id=submission_id)
        checksum_logger.set_state(file_path, file_metadata, state)
        if file_metadata.file_type in ("fastq", "bam"):
            seq_data_logger.set_state(file_path, file_metadata, state)


def test_validation_skips_files_already_validated_for_same_submission(
    submission,
    temp_checksum_log,
    temp_seq_data_log,
    mocker,
):
    """All files passing validation for the current submission_id should skip grz-check."""
    _seed_validation_logs(
        submission,
        temp_checksum_log,
        temp_seq_data_log,
        submission_id=submission.submission_id,
        validation_passed=True,
    )

    mock_validate = mocker.patch("grz_common.workers.submission.grz_check")

    list(
        submission.validate_files(
            checksum_progress_file=temp_checksum_log,
            seq_data_progress_file=temp_seq_data_log,
            threads=1,
            no_mmap=True,
        )
    )

    assert not mock_validate.validate_fastq_paired.called
    assert not mock_validate.validate_fastq.called
    assert not mock_validate.validate_bam.called
    assert not mock_validate.validate_raw.called


def test_validation_reruns_for_different_submission_id(
    submission,
    temp_checksum_log,
    temp_seq_data_log,
    mocker,
):
    """Files validated for a different submission_id must be re-validated."""
    _seed_validation_logs(
        submission,
        temp_checksum_log,
        temp_seq_data_log,
        submission_id="different-submission-id-9999",
        validation_passed=True,
    )

    mock_validate = mocker.patch("grz_common.workers.submission.grz_check")
    mock_report = mocker.MagicMock()
    mock_report.warnings = []
    mock_report.errors = []
    mock_report.is_valid = True
    mock_report.sha256 = None
    mock_validate.validate_fastq_paired.return_value = [mock_report, mock_report]
    mock_validate.validate_fastq.return_value = mock_report
    mock_validate.validate_bam.return_value = mock_report
    mock_validate.validate_raw.return_value = mock_report

    list(
        submission.validate_files(
            checksum_progress_file=temp_checksum_log,
            seq_data_progress_file=temp_seq_data_log,
            threads=1,
            no_mmap=True,
        )
    )

    assert (
        mock_validate.validate_fastq_paired.called
        or mock_validate.validate_fastq.called
        or mock_validate.validate_bam.called
        or mock_validate.validate_raw.called
    ), "Expected grz-check to re-run for a different submission_id"


def test_validation_reruns_after_failed_validation(
    submission,
    temp_checksum_log,
    temp_seq_data_log,
    mocker,
):
    """Files with validation_passed=False must be re-validated even with matching submission_id."""
    _seed_validation_logs(
        submission,
        temp_checksum_log,
        temp_seq_data_log,
        submission_id=submission.submission_id,
        validation_passed=False,
    )

    mock_validate = mocker.patch("grz_common.workers.submission.grz_check")
    mock_report = mocker.MagicMock()
    mock_report.warnings = []
    mock_report.errors = []
    mock_report.is_valid = True
    mock_report.sha256 = None
    mock_validate.validate_fastq_paired.return_value = [mock_report, mock_report]
    mock_validate.validate_fastq.return_value = mock_report
    mock_validate.validate_bam.return_value = mock_report
    mock_validate.validate_raw.return_value = mock_report

    list(
        submission.validate_files(
            checksum_progress_file=temp_checksum_log,
            seq_data_progress_file=temp_seq_data_log,
            threads=1,
            no_mmap=True,
        )
    )

    assert (
        mock_validate.validate_fastq_paired.called
        or mock_validate.validate_fastq.called
        or mock_validate.validate_bam.called
        or mock_validate.validate_raw.called
    ), "Expected grz-check to re-run for previously failed files"


def test_an_interruption_cancels_the_queued_validation_tasks(
    submission,
    temp_checksum_log,
    temp_seq_data_log,
    mocker,
):
    """A KeyboardInterrupt stops the validation without running the queued tasks. grzctl raises one for SIGTERM."""
    started_tasks = []

    def interrupt_the_first_task(*_args, **_kwargs):
        started_tasks.append(None)
        if len(started_tasks) == 1:
            # future.result() raises it again in the main thread, where a signal would raise it
            raise KeyboardInterrupt
        time.sleep(0.1)

    mock_validate = mocker.patch("grz_common.workers.submission.grz_check")
    mock_validate.validate_fastq_paired.side_effect = interrupt_the_first_task
    mock_validate.validate_fastq.side_effect = interrupt_the_first_task
    mock_validate.validate_bam.side_effect = interrupt_the_first_task
    mock_validate.validate_raw.side_effect = interrupt_the_first_task

    with pytest.raises(KeyboardInterrupt):
        list(
            submission.validate_files(
                checksum_progress_file=temp_checksum_log,
                seq_data_progress_file=temp_seq_data_log,
                threads=1,
                no_mmap=True,
            )
        )

    # the one worker thread may have started the second task before the interruption reached the main thread
    assert len(started_tasks) <= 2, f"{len(started_tasks)} tasks ran"
