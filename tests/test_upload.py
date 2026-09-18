"""Tests for the upload module"""

from pathlib import Path

import pytest
from grz_common.constants import MULTIPART_MIN_PART_SIZE
from grz_common.progress.progress_logging import FileProgressLogger
from grz_common.progress.states import UploadState
from grz_common.utils.checksums import calculate_sha256
from grz_common.workers.upload import S3BotoUploadWorker


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


def test_boto_upload(
    s3_config_model,
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
        s3_options=s3_config_model.s3,
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


def test_boto_upload_uses_configured_part_size(s3_config_model, remote_bucket, temp_upload_log_file_path, tmp_path):
    """upload_file splits a file into parts of the configured multipart_chunksize."""
    s3_options = s3_config_model.s3.model_copy(update={"multipart_chunksize": MULTIPART_MIN_PART_SIZE})
    upload_worker = S3BotoUploadWorker(s3_options=s3_options, status_file_path=temp_upload_log_file_path)
    file_path = tmp_path / "three_parts.bin"
    file_path.write_bytes(b"x" * (2 * MULTIPART_MIN_PART_SIZE + 1))

    upload_worker.upload_file(file_path, "three_parts.bin")

    etag = remote_bucket.Object("three_parts.bin").e_tag.strip('"')
    assert etag.endswith("-3"), f"expected a 3-part multipart upload, got ETag {etag}"


def test__gather_files_to_upload(encrypted_submission):
    submission_id = encrypted_submission.submission_id
    metadata_file_path, metadata_s3_object_id = encrypted_submission.get_metadata_file_path_and_object_id()
    gathered_files = encrypted_submission.get_encrypted_files_and_object_id()
    gathered_files[metadata_file_path] = metadata_s3_object_id
    gathered_files = sorted([(str(key), str(value)) for key, value in gathered_files.items()])

    expected_files = [
        (
            "tests/mock_files/submissions/valid_submission/encrypted_files/target_regions.bed.c4gh",
            f"{submission_id}/files/target_regions.bed.c4gh",
        ),
        (
            "tests/mock_files/submissions/valid_submission/encrypted_files/aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000_blood_normal.read1.fastq.gz.c4gh",
            f"{submission_id}/files/aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000_blood_normal.read1.fastq.gz.c4gh",
        ),
        (
            "tests/mock_files/submissions/valid_submission/encrypted_files/aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000_blood_normal.read2.fastq.gz.c4gh",
            f"{submission_id}/files/aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000_blood_normal.read2.fastq.gz.c4gh",
        ),
        (
            "tests/mock_files/submissions/valid_submission/encrypted_files/aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000_blood_normal.vcf.c4gh",
            f"{submission_id}/files/aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000_blood_normal.vcf.c4gh",
        ),
        (
            "tests/mock_files/submissions/valid_submission/encrypted_files/aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000_blood_tumor.read1.fastq.gz.c4gh",
            f"{submission_id}/files/aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000_blood_tumor.read1.fastq.gz.c4gh",
        ),
        (
            "tests/mock_files/submissions/valid_submission/encrypted_files/aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000_blood_tumor.read2.fastq.gz.c4gh",
            f"{submission_id}/files/aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000_blood_tumor.read2.fastq.gz.c4gh",
        ),
        (
            "tests/mock_files/submissions/valid_submission/encrypted_files/aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000_blood_tumor.vcf.c4gh",
            f"{submission_id}/files/aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000_blood_tumor.vcf.c4gh",
        ),
        (
            "tests/mock_files/submissions/valid_submission/encrypted_files/bbbbbbbb11111111bbbbbbbb11111111bbbbbbbb11111111bbbbbbbb11111111_blood_normal.read1.fastq.gz.c4gh",
            f"{submission_id}/files/bbbbbbbb11111111bbbbbbbb11111111bbbbbbbb11111111bbbbbbbb11111111_blood_normal.read1.fastq.gz.c4gh",
        ),
        (
            "tests/mock_files/submissions/valid_submission/encrypted_files/bbbbbbbb11111111bbbbbbbb11111111bbbbbbbb11111111bbbbbbbb11111111_blood_normal.read2.fastq.gz.c4gh",
            f"{submission_id}/files/bbbbbbbb11111111bbbbbbbb11111111bbbbbbbb11111111bbbbbbbb11111111_blood_normal.read2.fastq.gz.c4gh",
        ),
        (
            "tests/mock_files/submissions/valid_submission/encrypted_files/bbbbbbbb11111111bbbbbbbb11111111bbbbbbbb11111111bbbbbbbb11111111_blood_normal.vcf.c4gh",
            f"{submission_id}/files/bbbbbbbb11111111bbbbbbbb11111111bbbbbbbb11111111bbbbbbbb11111111_blood_normal.vcf.c4gh",
        ),
        (
            "tests/mock_files/submissions/valid_submission/metadata/metadata.json",
            f"{submission_id}/metadata/metadata.json",
        ),
    ]
    expected_files = sorted(expected_files)
    assert gathered_files == expected_files


def test_upload_skips_file_already_uploaded_for_same_submission(
    s3_config_model,
    remote_bucket,
    encrypted_submission,
    temp_upload_log_file_path,
    mocker,
):
    """Files with upload_successful=True and matching submission_id should be skipped."""
    upload_worker = S3BotoUploadWorker(
        s3_options=s3_config_model.s3,
        status_file_path=temp_upload_log_file_path,
    )

    progress_logger = FileProgressLogger[UploadState](temp_upload_log_file_path)
    for file_path, file_metadata in encrypted_submission.encrypted_files.items():
        progress_logger.set_state(
            file_path,
            file_metadata,
            state=UploadState(upload_successful=True, submission_id=encrypted_submission.submission_id),
        )

    files_to_upload = encrypted_submission.get_encrypted_files_and_object_id()
    metadata_path, metadata_s3_object_id = encrypted_submission.get_metadata_file_path_and_object_id()
    files_to_upload[metadata_path] = metadata_s3_object_id

    upload_spy = mocker.spy(upload_worker, "upload_file")
    upload_worker._upload_logged_files(encrypted_submission, progress_logger, files_to_upload)

    assert upload_spy.call_count == 0, f"Expected all files to be skipped, but {upload_spy.call_count} were uploaded"


def test_upload_rereuploads_file_with_different_submission_id(
    s3_config_model,
    remote_bucket,
    encrypted_submission,
    temp_upload_log_file_path,
    mocker,
):
    """Files logged as upload_successful=True for a different submission_id must be re-uploaded."""
    upload_worker = S3BotoUploadWorker(
        s3_options=s3_config_model.s3,
        status_file_path=temp_upload_log_file_path,
    )

    progress_logger = FileProgressLogger[UploadState](temp_upload_log_file_path)
    for file_path, file_metadata in encrypted_submission.encrypted_files.items():
        progress_logger.set_state(
            file_path,
            file_metadata,
            state=UploadState(upload_successful=True, submission_id="different-submission-id-9999"),
        )

    files_to_upload = encrypted_submission.get_encrypted_files_and_object_id()
    metadata_path, metadata_s3_object_id = encrypted_submission.get_metadata_file_path_and_object_id()
    files_to_upload[metadata_path] = metadata_s3_object_id

    upload_spy = mocker.spy(upload_worker, "upload_file")
    upload_worker._upload_logged_files(encrypted_submission, progress_logger, files_to_upload)

    expected = len(encrypted_submission.encrypted_files)
    assert upload_spy.call_count == expected, (
        f"Expected {expected} files to be re-uploaded for a different submission_id, "
        f"but only {upload_spy.call_count} were uploaded"
    )


def test_upload_rereuploads_file_after_failed_upload(
    s3_config_model,
    remote_bucket,
    encrypted_submission,
    temp_upload_log_file_path,
    mocker,
):
    """Files logged as upload_successful=False must be retried even with matching submission_id."""
    upload_worker = S3BotoUploadWorker(
        s3_options=s3_config_model.s3,
        status_file_path=temp_upload_log_file_path,
    )

    progress_logger = FileProgressLogger[UploadState](temp_upload_log_file_path)
    for file_path, file_metadata in encrypted_submission.encrypted_files.items():
        progress_logger.set_state(
            file_path,
            file_metadata,
            state=UploadState(upload_successful=False, submission_id=encrypted_submission.submission_id),
        )

    files_to_upload = encrypted_submission.get_encrypted_files_and_object_id()
    metadata_path, metadata_s3_object_id = encrypted_submission.get_metadata_file_path_and_object_id()
    files_to_upload[metadata_path] = metadata_s3_object_id

    upload_spy = mocker.spy(upload_worker, "upload_file")
    upload_worker._upload_logged_files(encrypted_submission, progress_logger, files_to_upload)

    expected = len(encrypted_submission.encrypted_files)
    assert upload_spy.call_count == expected, (
        f"Expected {expected} failed files to be retried, but only {upload_spy.call_count} were uploaded"
    )


def test_upload_file_uploads_parts_in_parallel(
    s3_config_model, remote_bucket, temp_upload_log_file_path, overlapping_s3_calls, tmp_path
):
    """With two threads, two parts of one file are uploaded at the same time."""
    s3_options = s3_config_model.s3.model_copy(update={"multipart_chunksize": MULTIPART_MIN_PART_SIZE})
    file_path = tmp_path / "three_parts.bin"
    file_path.write_bytes(b"x" * (2 * MULTIPART_MIN_PART_SIZE + 1))
    overlapped = overlapping_s3_calls("UploadPart", timeout=10)
    upload_worker = S3BotoUploadWorker(s3_options=s3_options, status_file_path=temp_upload_log_file_path, threads=2)

    upload_worker.upload_file(file_path, "three_parts.bin")

    assert overlapped(), "two parts should have been uploaded at the same time"


def test_upload_file_uploads_one_part_at_a_time_with_one_thread(
    s3_config_model, remote_bucket, temp_upload_log_file_path, overlapping_s3_calls, tmp_path
):
    """With one thread, no two parts of a file are uploaded at the same time."""
    s3_options = s3_config_model.s3.model_copy(update={"multipart_chunksize": MULTIPART_MIN_PART_SIZE})
    file_path = tmp_path / "three_parts.bin"
    file_path.write_bytes(b"x" * (2 * MULTIPART_MIN_PART_SIZE + 1))
    overlapped = overlapping_s3_calls("UploadPart", timeout=0.5)
    upload_worker = S3BotoUploadWorker(s3_options=s3_options, status_file_path=temp_upload_log_file_path, threads=1)

    upload_worker.upload_file(file_path, "three_parts.bin")

    assert not overlapped(), "one thread should upload one part after the other"
