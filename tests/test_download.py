"""Tests for the download module"""

import io
from pathlib import Path
from unittest.mock import MagicMock

import botocore.client
import pytest
from grz_common.exceptions import DownloadError, NetworkError
from grz_common.utils.checksums import calculate_sha256
from grz_common.workers.download import S3BotoDownloadWorker
from grz_common.workers.submission import EncryptedSubmission
from grz_common.workers.worker import Worker


class _BreakingBody(io.RawIOBase):
    """An S3 response body that hands out one chunk and then breaks, like a dropped connection."""

    def __init__(self, first_chunk: bytes):
        super().__init__()
        self._first_chunk = first_chunk

    def readable(self) -> bool:
        return True

    def read(self, size: int | None = -1) -> bytes:
        if not self._first_chunk:
            raise ConnectionError("connection reset by peer")
        chunk = self._first_chunk
        self._first_chunk = b""
        return chunk


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

    mock_meta = MagicMock()

    # Execute download
    local_file_path = files_dir / "large_test_file.fastq"
    mock_meta.file_size_in_bytes = 1000
    s3_object_id = f"{submission_id}/large_test_file.fastq"
    download_worker._download_with_progress(str(local_file_path), s3_object_id, mock_meta)

    local_file_path = files_dir / "small_test_file.txt"
    mock_meta.file_size_in_bytes = 100
    s3_object_id = f"{submission_id}/small_test_file.txt"
    download_worker._download_with_progress(str(local_file_path), s3_object_id, mock_meta)

    # Assert that the files have been downloaded correctly
    assert (files_dir / "large_test_file.fastq").exists(), "Fastq file was not downloaded."
    assert (files_dir / "small_test_file.txt").exists(), "Text file was not downloaded."

    # Further assertions can be made here as necessary
    assert calculate_sha256(files_dir / "large_test_file.fastq") == temp_fastq_file_sha256sum, "Fastq SHA256 mismatch."
    assert calculate_sha256(files_dir / "small_test_file.txt") == temp_small_file_sha256sum, (
        "Text file SHA256 mismatch."
    )


def test_download_file_fails_for_missing_key(
    s3_config_model,
    remote_bucket,
    encrypted_submission,
    tmp_path,
):
    """A key that is not in the bucket fails with a DownloadError."""
    from grz_common.progress.progress_logging import FileProgressLogger
    from grz_common.progress.states import DownloadState

    download_log_path = tmp_path / "progress_download.cjson"
    download_worker = S3BotoDownloadWorker(
        s3_options=s3_config_model.s3,
        status_file_path=download_log_path,
    )
    progress_logger = FileProgressLogger[DownloadState](download_log_path)
    file_path, file_metadata = next(iter(encrypted_submission.encrypted_files.items()))

    with pytest.raises(DownloadError):
        download_worker.download_file(
            tmp_path / "files" / file_path.name,
            f"{encrypted_submission.submission_id}/files/missing.c4gh",
            progress_logger,
            file_metadata,
            encrypted_submission.submission_id,
        )


def test_download_file_leaves_no_file_behind_for_a_missing_key(
    s3_config_model,
    remote_bucket,
    encrypted_submission,
    tmp_path,
):
    """A key that is not in the bucket leaves no local file behind."""
    from grz_common.progress.progress_logging import FileProgressLogger
    from grz_common.progress.states import DownloadState

    download_log_path = tmp_path / "progress_download.cjson"
    download_worker = S3BotoDownloadWorker(
        s3_options=s3_config_model.s3,
        status_file_path=download_log_path,
    )
    progress_logger = FileProgressLogger[DownloadState](download_log_path)
    file_path, file_metadata = next(iter(encrypted_submission.encrypted_files.items()))
    local_file_path = tmp_path / "files" / file_path.name

    with pytest.raises(DownloadError):
        download_worker.download_file(
            local_file_path,
            f"{encrypted_submission.submission_id}/files/missing.c4gh",
            progress_logger,
            file_metadata,
            encrypted_submission.submission_id,
        )

    assert not local_file_path.exists(), "a failed download should leave no local file behind"


def test_download_file_leaves_no_partial_file_when_the_stream_breaks(
    s3_config_model,
    remote_bucket,
    submission_metadata_dir,
    monkeypatch,
    tmp_path,
):
    """A download that breaks part way through leaves no partial file behind."""
    from grz_common.progress.progress_logging import FileProgressLogger
    from grz_common.progress.states import DownloadState

    submission = _submission_in_the_bucket(remote_bucket, submission_metadata_dir, tmp_path)
    original_call = botocore.client.BaseClient._make_api_call

    def break_the_body(self, operation_name, kwargs):
        response = original_call(self, operation_name, kwargs)
        if operation_name == "GetObject":
            response["Body"] = _BreakingBody(b"encrypted ")
        return response

    monkeypatch.setattr(botocore.client.BaseClient, "_make_api_call", break_the_body)

    download_log_path = tmp_path / "progress_download.cjson"
    download_worker = S3BotoDownloadWorker(
        s3_options=s3_config_model.s3,
        status_file_path=download_log_path,
    )
    progress_logger = FileProgressLogger[DownloadState](download_log_path)
    local_file_path, file_metadata = next(iter(submission.encrypted_files.items()))

    with pytest.raises(NetworkError):
        download_worker.download_file(
            local_file_path,
            f"{submission.submission_id}/files/{file_metadata.encrypted_file_path()}",
            progress_logger,
            file_metadata,
            submission.submission_id,
        )

    assert not local_file_path.exists(), "a broken download should leave no partial file behind"


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


def _submission_in_the_bucket(remote_bucket, submission_metadata_dir: Path, tmp_path: Path) -> EncryptedSubmission:
    """Put every encrypted file of the example submission into the bucket, to be downloaded to *tmp_path*."""
    submission = EncryptedSubmission(submission_metadata_dir, tmp_path / "encrypted_files")
    for file_metadata in submission.encrypted_files.values():
        key = f"{submission.submission_id}/files/{file_metadata.encrypted_file_path()}"
        remote_bucket.put_object(Key=key, Body=b"encrypted payload")
    return submission


def test_download_downloads_files_in_parallel(
    s3_config_model, remote_bucket, submission_metadata_dir, temp_download_log_file_path, overlapping_s3_calls, tmp_path
):
    """With two threads, two files of a submission are downloaded at the same time."""
    submission = _submission_in_the_bucket(remote_bucket, submission_metadata_dir, tmp_path)
    overlapped = overlapping_s3_calls("GetObject", timeout=10)
    download_worker = S3BotoDownloadWorker(
        s3_options=s3_config_model.s3, status_file_path=temp_download_log_file_path, threads=2
    )

    download_worker.download(submission.submission_id, submission)

    assert overlapped(), "two files should have been downloaded at the same time"
    assert all(path.exists() for path in submission.encrypted_files)


def test_download_downloads_one_file_at_a_time_with_one_thread(
    s3_config_model, remote_bucket, submission_metadata_dir, temp_download_log_file_path, overlapping_s3_calls, tmp_path
):
    """With one thread, no two files of a submission are downloaded at the same time."""
    submission = _submission_in_the_bucket(remote_bucket, submission_metadata_dir, tmp_path)
    overlapped = overlapping_s3_calls("GetObject", timeout=0.5)
    download_worker = S3BotoDownloadWorker(
        s3_options=s3_config_model.s3, status_file_path=temp_download_log_file_path, threads=1
    )

    download_worker.download(submission.submission_id, submission)

    assert not overlapped(), "one thread should download one file after the other"
    assert all(path.exists() for path in submission.encrypted_files)
