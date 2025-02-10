"""Tests for the download module"""

from pathlib import Path

import pytest
from moto import mock_aws

from grz_cli.download import S3BotoDownloadWorker
from grz_cli.file_operations import calculate_sha256


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


@mock_aws
def test_boto_download(
    config_model,
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
        config=config_model,
        status_file_path=temp_download_log_file_path,
    )

    # Execute download
    local_file_path = files_dir / "large_test_file.fastq"
    s3_object_id = f"{submission_id}/large_test_file.fastq"
    download_worker.download_file(local_file_path, s3_object_id, 100000)
    local_file_path = files_dir / "small_test_file.txt"
    s3_object_id = f"{submission_id}/small_test_file.txt"
    download_worker.download_file(local_file_path, s3_object_id, 100000)

    # Assert that the files have been downloaded correctly
    assert (files_dir / "large_test_file.fastq").exists(), "Fastq file was not downloaded."
    assert (files_dir / "small_test_file.txt").exists(), "Text file was not downloaded."

    # Further assertions can be made here as necessary
    assert calculate_sha256(files_dir / "large_test_file.fastq") == temp_fastq_file_sha256sum, "Fastq SHA256 mismatch."
    assert calculate_sha256(files_dir / "small_test_file.txt") == temp_small_file_sha256sum, (
        "Text file SHA256 mismatch."
    )
