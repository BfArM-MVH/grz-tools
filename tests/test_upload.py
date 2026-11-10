"""Tests for the upload module"""

from http import HTTPStatus
from pathlib import Path

import botocore.client
import grz_common.exceptions as grzexc
import pytest
from botocore.exceptions import ClientError
from grz_common.progress.progress_logging import FileProgressLogger
from grz_common.progress.states import UploadState
from grz_common.utils.checksums import calculate_sha256
from grz_common.workers.submission import EncryptedSubmission
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


def test_upload_of_a_submission_already_in_the_inbox_raises(
    s3_config_model, remote_bucket, encrypted_submission, tmp_path
):
    """An inbox that already holds the submission's metadata means that the submission was uploaded before."""
    _, metadata_s3_object_id = encrypted_submission.get_metadata_file_path_and_object_id()
    remote_bucket.put_object(Key=metadata_s3_object_id, Body=b"{}")
    upload_worker = S3BotoUploadWorker(
        s3_options=s3_config_model.s3, status_file_path=tmp_path / "progress_upload.cjson"
    )

    with pytest.raises(grzexc.DuplicateUploadError):
        upload_worker.upload(encrypted_submission)


def test_archive_of_an_archived_submission_uploads_nothing(
    s3_config_model, remote_bucket, encrypted_submission, tmp_path, mocker
):
    """Archive uploads the metadata last, so a rerun after a finished archival returns without an upload."""
    upload_worker = S3BotoUploadWorker(
        s3_options=s3_config_model.s3, status_file_path=tmp_path / "progress_upload.cjson"
    )
    upload_worker.archive(encrypted_submission)
    _, metadata_s3_object_id = encrypted_submission.get_metadata_file_path_and_object_id()
    assert metadata_s3_object_id in {o.key for o in remote_bucket.objects.all()}

    upload_spy = mocker.spy(upload_worker, "upload_file")
    upload_worker.archive(encrypted_submission)

    assert upload_spy.call_count == 0, "a rerun must not upload anything"


def test_upload_of_a_submission_with_a_missing_local_file_raises(
    s3_config_model, remote_bucket, submission_metadata_dir, tmp_path
):
    """A file that the encryption did not produce is a step that has not passed, not a failed upload."""
    submission = EncryptedSubmission(submission_metadata_dir, tmp_path / "encrypted_files")
    upload_worker = S3BotoUploadWorker(
        s3_options=s3_config_model.s3, status_file_path=tmp_path / "progress_upload.cjson"
    )

    with pytest.raises(grzexc.IncompleteSubmissionError):
        upload_worker.upload(submission)


def _fail_s3_operation(monkeypatch, operation: str, code: str, status: HTTPStatus = HTTPStatus.FORBIDDEN):
    """Answer every S3 call of *operation* with the error *code* and the HTTP *status*."""
    original_call = botocore.client.BaseClient._make_api_call

    def fail(self, operation_name, kwargs):
        if operation_name == operation:
            raise ClientError(
                {"Error": {"Code": code, "Message": code}, "ResponseMetadata": {"HTTPStatusCode": status.value}},
                operation_name,
            )
        return original_call(self, operation_name, kwargs)

    monkeypatch.setattr(botocore.client.BaseClient, "_make_api_call", fail)


def test_upload_file_reports_rejected_credentials_as_a_configuration_error(
    s3_config_model, remote_bucket, temp_small_file_path, tmp_path, monkeypatch
):
    """S3Transfer wraps the ClientError in S3UploadFailedError, and the upload worker takes it out again."""
    _fail_s3_operation(monkeypatch, "PutObject", "InvalidAccessKeyId")
    upload_worker = S3BotoUploadWorker(
        s3_options=s3_config_model.s3, status_file_path=tmp_path / "progress_upload.cjson"
    )

    with pytest.raises(grzexc.ConfigurationError) as excinfo:
        upload_worker.upload_file(temp_small_file_path, "small_test_file.bed")

    assert isinstance(excinfo.value.__cause__, ClientError)


def test_upload_file_reports_any_other_s3_error_as_a_failed_upload(
    s3_config_model, remote_bucket, temp_small_file_path, tmp_path, monkeypatch
):
    """``AccessDenied`` is no configuration error, since S3 also answers it for a missing object."""
    _fail_s3_operation(monkeypatch, "PutObject", "AccessDenied")
    upload_worker = S3BotoUploadWorker(
        s3_options=s3_config_model.s3, status_file_path=tmp_path / "progress_upload.cjson"
    )

    with pytest.raises(grzexc.UploadError):
        upload_worker.upload_file(temp_small_file_path, "small_test_file.bed")


def test_upload_reports_a_faulty_setup_as_a_configuration_error(
    s3_config_model, faulty_s3_setup, encrypted_submission, tmp_path
):
    """The check for an earlier upload is the first request, and it already finds a faulty setup."""
    upload_worker = S3BotoUploadWorker(
        s3_options=s3_config_model.s3, status_file_path=tmp_path / "progress_upload.cjson"
    )

    with pytest.raises(grzexc.ConfigurationError):
        upload_worker.upload(encrypted_submission)


@pytest.mark.parametrize("method", ["upload", "archive"])
def test_upload_goes_ahead_with_a_warning_if_access_is_denied(
    s3_config_model, remote_bucket, encrypted_submission, tmp_path, monkeypatch, caplog, method
):
    """Credentials that may not list the bucket get ``AccessDenied`` for a missing object, so the check cannot tell."""
    _fail_s3_operation(monkeypatch, "HeadObject", "403")
    _fail_s3_operation(monkeypatch, "GetObject", "AccessDenied")
    upload_worker = S3BotoUploadWorker(
        s3_options=s3_config_model.s3, status_file_path=tmp_path / "progress_upload.cjson"
    )

    getattr(upload_worker, method)(encrypted_submission)

    _, metadata_s3_object_id = encrypted_submission.get_metadata_file_path_and_object_id()
    assert metadata_s3_object_id in {o.key for o in remote_bucket.objects.all()}
    assert "because S3 denies access" in caplog.text
