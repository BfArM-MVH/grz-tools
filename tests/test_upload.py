"""Tests for the upload module"""

from pathlib import Path

import pytest
from moto import mock_aws

from grz_cli.file_operations import calculate_sha256
from grz_cli.upload import S3BotoUploadWorker


@pytest.fixture(scope="module")
def temp_log_dir(tmpdir_factory: pytest.TempdirFactory):
    """Create temporary log folder for this pytest module"""
    datadir = tmpdir_factory.mktemp("logs")
    return datadir


@pytest.fixture
def temp_upload_log_file_path(temp_log_dir) -> Path:
    log_file = Path(temp_log_dir) / "progress_upload.cjson"
    return log_file


def download_file(remote_bucket, object_id, output_path):
    remote_bucket.download_file(object_id, output_path)


@mock_aws
def test_boto_upload(
    config_model,
    remote_bucket,
    temp_small_file_path,
    temp_small_file_sha256sum,
    temp_fastq_file_path,
    temp_fastq_file_sha256sum,
    temp_upload_log_file_path,
    tmpdir_factory,
):
    # create upload worker
    upload_worker = S3BotoUploadWorker(
        config=config_model,
        status_file_path=temp_upload_log_file_path,
    )

    upload_worker.upload_file(temp_small_file_path, "small_test_file.bed")
    upload_worker.upload_file(temp_fastq_file_path, "large_test_file.fastq")

    # download files again
    local_tmpdir = tmpdir_factory.mktemp("redownload")
    local_tmpdir_path = Path(local_tmpdir.strpath)

    download_file(remote_bucket, "small_test_file.bed", local_tmpdir_path / "small_test_file.bed")
    download_file(
        remote_bucket,
        "large_test_file.fastq",
        local_tmpdir_path / "large_test_file.fastq",
    )

    assert calculate_sha256(local_tmpdir_path / "small_test_file.bed") == temp_small_file_sha256sum
    assert calculate_sha256(local_tmpdir_path / "large_test_file.fastq") == temp_fastq_file_sha256sum


def test__gather_files_to_upload(encrypted_submission):
    metadata_file_path, metadata_s3_object_id = encrypted_submission.get_metadata_file_path_and_object_id()
    gathered_files = encrypted_submission.get_encrypted_files_and_object_id()
    gathered_files[metadata_file_path] = metadata_s3_object_id
    gathered_files = sorted([(str(key), str(value)) for key, value in gathered_files.items()])

    expected_files = [
        (
            "tests/mock_files/submissions/valid_submission/encrypted_files/target_regions.bed.c4gh",
            "aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000/files/target_regions.bed.c4gh",
        ),
        (
            "tests/mock_files/submissions/valid_submission/encrypted_files/aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000_blood_normal.read1.fastq.gz.c4gh",
            "aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000/files/aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000_blood_normal.read1.fastq.gz.c4gh",
        ),
        (
            "tests/mock_files/submissions/valid_submission/encrypted_files/aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000_blood_normal.read2.fastq.gz.c4gh",
            "aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000/files/aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000_blood_normal.read2.fastq.gz.c4gh",
        ),
        (
            "tests/mock_files/submissions/valid_submission/encrypted_files/aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000_blood_normal.vcf.c4gh",
            "aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000/files/aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000_blood_normal.vcf.c4gh",
        ),
        (
            "tests/mock_files/submissions/valid_submission/encrypted_files/aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000_blood_tumor.read1.fastq.gz.c4gh",
            "aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000/files/aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000_blood_tumor.read1.fastq.gz.c4gh",
        ),
        (
            "tests/mock_files/submissions/valid_submission/encrypted_files/aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000_blood_tumor.read2.fastq.gz.c4gh",
            "aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000/files/aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000_blood_tumor.read2.fastq.gz.c4gh",
        ),
        (
            "tests/mock_files/submissions/valid_submission/encrypted_files/aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000_blood_tumor.vcf.c4gh",
            "aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000/files/aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000_blood_tumor.vcf.c4gh",
        ),
        (
            "tests/mock_files/submissions/valid_submission/encrypted_files/bbbbbbbb11111111bbbbbbbb11111111bbbbbbbb11111111bbbbbbbb11111111_blood_normal.read1.fastq.gz.c4gh",
            "aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000/files/bbbbbbbb11111111bbbbbbbb11111111bbbbbbbb11111111bbbbbbbb11111111_blood_normal.read1.fastq.gz.c4gh",
        ),
        (
            "tests/mock_files/submissions/valid_submission/encrypted_files/bbbbbbbb11111111bbbbbbbb11111111bbbbbbbb11111111bbbbbbbb11111111_blood_normal.read2.fastq.gz.c4gh",
            "aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000/files/bbbbbbbb11111111bbbbbbbb11111111bbbbbbbb11111111bbbbbbbb11111111_blood_normal.read2.fastq.gz.c4gh",
        ),
        (
            "tests/mock_files/submissions/valid_submission/encrypted_files/bbbbbbbb11111111bbbbbbbb11111111bbbbbbbb11111111bbbbbbbb11111111_blood_normal.vcf.c4gh",
            "aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000/files/bbbbbbbb11111111bbbbbbbb11111111bbbbbbbb11111111bbbbbbbb11111111_blood_normal.vcf.c4gh",
        ),
        (
            "tests/mock_files/submissions/valid_submission/metadata/metadata.json",
            "aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000/metadata/metadata.json",
        ),
    ]
    expected_files = sorted(expected_files)
    assert gathered_files == expected_files
