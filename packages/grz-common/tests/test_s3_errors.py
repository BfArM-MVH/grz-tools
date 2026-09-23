"""An error of the S3 client becomes the failure it stands for."""

import grz_common.exceptions as grzexc
import pytest
from boto3.exceptions import S3UploadFailedError
from botocore.exceptions import ClientError, EndpointConnectionError, NoCredentialsError
from grz_common.transfer import get_metadata_upload_timestamp, head_object, s3_error


def _client_error(code: str) -> ClientError:
    return ClientError({"Error": {"Code": code, "Message": code}}, "HeadObject")


class _FailingClient:
    def __init__(self, error: Exception):
        self.error = error

    def head_object(self, **_kwargs):
        raise self.error


@pytest.mark.parametrize("code", ["InvalidAccessKeyId", "SignatureDoesNotMatch", "NoSuchBucket"])
def test_an_error_that_only_a_faulty_setup_causes_is_a_configuration_error(code: str):
    assert isinstance(s3_error(_client_error(code), "Reading s3://bucket/key"), grzexc.ConfigurationError)


def test_missing_credentials_are_a_configuration_error():
    assert isinstance(s3_error(NoCredentialsError(), "Reading s3://bucket/key"), grzexc.ConfigurationError)


def _upload_failed_while_handling(error: ClientError) -> S3UploadFailedError:
    """Raise ``S3UploadFailedError`` the way ``S3Transfer.upload_file`` does: inside ``except``, without ``from``."""
    try:
        try:
            raise error
        except ClientError as e:
            raise S3UploadFailedError(f"Failed to upload file to bucket/key: {e}")  # noqa: B904
    except S3UploadFailedError as wrapped:
        return wrapped


def test_the_error_code_counts_inside_the_wrapper_of_s3transfer():
    failure = s3_error(_upload_failed_while_handling(_client_error("InvalidAccessKeyId")), "Upload", grzexc.UploadError)

    assert isinstance(failure, grzexc.ConfigurationError)


def test_the_wrapper_of_s3transfer_is_a_failed_transfer_for_any_other_code():
    failure = s3_error(_upload_failed_while_handling(_client_error("SlowDown")), "Upload", grzexc.UploadError)

    assert type(failure) is grzexc.UploadError


def test_a_suppressed_context_does_not_count():
    """``raise ... from None`` hides the context, so its error code does not count either."""
    try:
        try:
            raise _client_error("InvalidAccessKeyId")
        except ClientError:
            raise S3UploadFailedError("upload failed") from None
    except S3UploadFailedError as e:
        failure = s3_error(e, "Upload", grzexc.UploadError)

    assert type(failure) is grzexc.UploadError


@pytest.mark.parametrize(
    "error",
    [
        # also the answer for a missing object if the credentials may not list the bucket
        _client_error("AccessDenied"),
        _client_error("SlowDown"),
        EndpointConnectionError(endpoint_url="https://s3.localhost"),
    ],
    ids=["access-denied", "slow-down", "unreachable"],
)
def test_any_other_error_is_a_failed_transfer_of_the_given_kind(error: Exception):
    failure = s3_error(error, "Upload to s3://bucket/key", grzexc.UploadError)

    assert type(failure) is grzexc.UploadError
    assert str(failure).startswith("Upload to s3://bucket/key failed: ")


def test_head_object_reports_a_missing_object():
    with pytest.raises(grzexc.MissingObjectError):
        head_object(_FailingClient(_client_error("404")), "bucket", "key")


def test_head_object_reports_a_refused_read_as_a_failed_download():
    with pytest.raises(grzexc.DownloadError) as excinfo:
        head_object(_FailingClient(_client_error("AccessDenied")), "bucket", "key")

    assert not isinstance(excinfo.value, grzexc.MissingObjectError)


def test_a_missing_metadata_object_is_a_missing_submission_file():
    with pytest.raises(grzexc.MissingSubmissionFileError):
        get_metadata_upload_timestamp(_FailingClient(_client_error("404")), "bucket", "submission")
