"""Tests for the download module"""

from pathlib import Path

import pytest
from grz_common.utils.checksums import calculate_sha256
from grz_common.workers.download import S3BotoDownloadWorker
from grz_common.workers.worker import Worker


@pytest.fixture(scope="module")
def temp_log_dir(tmpdir_factory: pytest.TempdirFactory):
    """Create temporary log folder for this pytest module"""
    datadir = tmpdir_factory.mktemp("logs")
    return datadir


@pytest.fixture
def temp_download_log_file_path(temp_log_dir) -> Path:
    log_file = Path(temp_log_dir) / "progress_download.cjson"
    return log_file


def upload_file(remote_bucket, local_file_path, s3_key):
    """Upload file to the remote S3 bucket with the specified key."""
    remote_bucket.upload_file(local_file_path, s3_key)


def test_boto_download(
    s3_config_model,
    remote_bucket,
    temp_small_file_path,
    temp_small_file_sha256sum,
    temp_fastq_file_path,
    temp_fastq_file_sha256sum,
    temp_download_log_file_path,
    tmpdir_factory,
):
    # Prepare directories
    submission_id = "submission123"  # Use the same submission ID as in the download method

    files_dir = Path(tmpdir_factory.mktemp(submission_id))

    # Upload metadata and file using the correct submission ID
    upload_file(remote_bucket, temp_fastq_file_path, f"{submission_id}/large_test_file.fastq")
    upload_file(remote_bucket, temp_small_file_path, f"{submission_id}/small_test_file.txt")

    # Create a mock S3 bucket
    download_worker = S3BotoDownloadWorker(
        s3_options=s3_config_model.s3,
        status_file_path=temp_download_log_file_path,
    )

    # Execute download
    local_file_path = files_dir / "large_test_file.fastq"
    s3_object_id = f"{submission_id}/large_test_file.fastq"
    download_worker._download_with_progress(str(local_file_path), s3_object_id)
    local_file_path = files_dir / "small_test_file.txt"
    s3_object_id = f"{submission_id}/small_test_file.txt"
    download_worker._download_with_progress(str(local_file_path), s3_object_id)

    # Assert that the files have been downloaded correctly
    assert (files_dir / "large_test_file.fastq").exists(), "Fastq file was not downloaded."
    assert (files_dir / "small_test_file.txt").exists(), "Text file was not downloaded."

    # Further assertions can be made here as necessary
    assert calculate_sha256(files_dir / "large_test_file.fastq") == temp_fastq_file_sha256sum, "Fastq SHA256 mismatch."
    assert calculate_sha256(files_dir / "small_test_file.txt") == temp_small_file_sha256sum, (
        "Text file SHA256 mismatch."
    )


def test_download_skips_file_already_downloaded_for_same_submission(
    s3_config_model,
    remote_bucket,
    temp_download_log_file_path,
    encrypted_submission,
    mocker,
):
    """Files with download_successful=True and matching submission_id should be skipped."""
    from grz_common.progress.progress_logging import FileProgressLogger
    from grz_common.progress.states import DownloadState

    download_worker = S3BotoDownloadWorker(
        s3_options=s3_config_model.s3,
        status_file_path=temp_download_log_file_path,
    )

    progress_logger = FileProgressLogger[DownloadState](temp_download_log_file_path)
    for file_path, file_metadata in encrypted_submission.encrypted_files.items():
        progress_logger.set_state(
            file_path,
            file_metadata,
            state=DownloadState(download_successful=True, submission_id=encrypted_submission.submission_id),
        )

    download_spy = mocker.spy(download_worker, "download_file")
    download_worker.download(encrypted_submission.submission_id, encrypted_submission)

    assert download_spy.call_count == 0, (
        f"Expected all files to be skipped, but {download_spy.call_count} were downloaded"
    )


def test_download_redownloads_file_with_different_submission_id(
    s3_config_model,
    remote_bucket,
    temp_download_log_file_path,
    encrypted_submission,
    mocker,
):
    """Files logged as download_successful=True for a different submission_id must be re-downloaded."""
    from grz_common.progress.progress_logging import FileProgressLogger
    from grz_common.progress.states import DownloadState

    download_worker = S3BotoDownloadWorker(
        s3_options=s3_config_model.s3,
        status_file_path=temp_download_log_file_path,
    )

    progress_logger = FileProgressLogger[DownloadState](temp_download_log_file_path)
    for file_path, file_metadata in encrypted_submission.encrypted_files.items():
        progress_logger.set_state(
            file_path,
            file_metadata,
            state=DownloadState(download_successful=True, submission_id="different-submission-id-9999"),
        )

    mock_download = mocker.patch.object(download_worker, "download_file")
    download_worker.download(encrypted_submission.submission_id, encrypted_submission)

    expected = len(encrypted_submission.encrypted_files)
    assert mock_download.call_count == expected, (
        f"Expected {expected} files to be re-downloaded for a different submission_id, "
        f"but only {mock_download.call_count} were downloaded"
    )


def test_download_redownloads_file_after_failed_download(
    s3_config_model,
    remote_bucket,
    temp_download_log_file_path,
    encrypted_submission,
    mocker,
):
    """Files logged as download_successful=False must be retried even with matching submission_id."""
    from grz_common.progress.progress_logging import FileProgressLogger
    from grz_common.progress.states import DownloadState

    download_worker = S3BotoDownloadWorker(
        s3_options=s3_config_model.s3,
        status_file_path=temp_download_log_file_path,
    )

    progress_logger = FileProgressLogger[DownloadState](temp_download_log_file_path)
    for file_path, file_metadata in encrypted_submission.encrypted_files.items():
        progress_logger.set_state(
            file_path,
            file_metadata,
            state=DownloadState(download_successful=False, submission_id=encrypted_submission.submission_id),
        )

    mock_download = mocker.patch.object(download_worker, "download_file")
    download_worker.download(encrypted_submission.submission_id, encrypted_submission)

    expected = len(encrypted_submission.encrypted_files)
    assert mock_download.call_count == expected, (
        f"Expected {expected} failed files to be retried, but only {mock_download.call_count} were downloaded"
    )


def test_worker_download_checks_metadata_version_before_files(
    s3_config_model,
    remote_bucket,
    encrypted_submission,
    tmp_path,
):
    """A failing metadata version check aborts before encrypted files are downloaded."""
    metadata_path, metadata_key = encrypted_submission.get_metadata_file_path_and_object_id()
    upload_file(remote_bucket, metadata_path, metadata_key)
    for local_file_path, s3_key in encrypted_submission.get_encrypted_files_and_object_id().items():
        upload_file(remote_bucket, local_file_path, s3_key)

    worker = Worker(
        metadata_dir=tmp_path / "metadata",
        files_dir=tmp_path / "files",
        log_dir=tmp_path / "logs",
        encrypted_files_dir=tmp_path / "encrypted_files",
    )

    checked_versions = []

    def fail_metadata_version_check(metadata_schema_version: str) -> None:
        checked_versions.append(metadata_schema_version)
        raise SystemExit(1)

    with pytest.raises(SystemExit):
        worker.download(
            s3_config_model.s3,
            encrypted_submission.submission_id,
            metadata_version_check=fail_metadata_version_check,
        )

    assert checked_versions == [encrypted_submission.metadata.content.get_schema_version()]
    assert (tmp_path / "metadata" / "metadata.json").exists()
    assert not list((tmp_path / "encrypted_files").rglob("*.c4gh"))
