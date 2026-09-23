"""An error of the S3 client becomes the failure it stands for."""

import io

import boto3
import grz_common.exceptions as grzexc
import pytest
from boto3.exceptions import S3UploadFailedError
from botocore.exceptions import ClientError, EndpointConnectionError, NoCredentialsError
from grz_common.transfer import get_metadata_upload_timestamp, head_object, s3_error


def _client_error(code: str, operation: str = "HeadObject") -> ClientError:
    return ClientError({"Error": {"Code": code, "Message": code}}, operation)


class _FailingClient:
    """Fail a HEAD request with the HTTP status alone, as S3 does, and a GET request with ``get_code``.

    Without ``get_code``, the GET request succeeds.
    """

    def __init__(self, head_status: str, get_code: str | None = None):
        self.head_status = head_status
        self.get_code = get_code
        self.body = io.BytesIO(b"x")

    def head_object(self, **_kwargs):
        raise _client_error(self.head_status)

    def get_object(self, **_kwargs):
        if self.get_code is None:
            return {"Body": self.body}
        raise _client_error(self.get_code, "GetObject")


@pytest.mark.parametrize("code", ["InvalidAccessKeyId", "SignatureDoesNotMatch", "NoSuchBucket"])
def test_an_error_that_only_a_faulty_setup_causes_is_a_configuration_error(code: str):
    assert isinstance(s3_error(_client_error(code), "Reading s3://bucket/key"), grzexc.ConfigurationError)


def test_missing_credentials_are_a_configuration_error():
    assert isinstance(s3_error(NoCredentialsError(), "Reading s3://bucket/key"), grzexc.ConfigurationError)


def test_a_bucket_name_that_botocore_rejects_is_a_configuration_error():
    """The S3 client checks the bucket name before it sends a request, so no S3 server is needed."""
    s3_client = boto3.client("s3", region_name="us-east-1", aws_access_key_id="key", aws_secret_access_key="secret")

    with pytest.raises(grzexc.ConfigurationError):
        head_object(s3_client, "not a bucket name", "key")


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
        head_object(_FailingClient("404", "NoSuchKey"), "bucket", "key")


@pytest.mark.parametrize(
    ("head_status", "get_code"),
    [("404", "NoSuchBucket"), ("403", "InvalidAccessKeyId"), ("403", "SignatureDoesNotMatch")],
)
def test_head_object_learns_a_faulty_setup_from_a_get_request(head_status: str, get_code: str):
    """S3 answers HEAD with the HTTP status alone, so the error code of a GET request decides."""
    with pytest.raises(grzexc.ConfigurationError):
        head_object(_FailingClient(head_status, get_code), "bucket", "key")


def test_head_object_reports_a_refused_read_as_a_failed_download():
    with pytest.raises(grzexc.DownloadError) as excinfo:
        head_object(_FailingClient("403", "AccessDenied"), "bucket", "key")

    assert not isinstance(excinfo.value, grzexc.MissingObjectError)


def test_head_object_raises_the_given_class_for_a_failed_transfer():
    with pytest.raises(grzexc.UploadError):
        head_object(_FailingClient("403", "AccessDenied"), "bucket", "key", grzexc.UploadError)


@pytest.mark.parametrize(
    ("head_status", "expected"),
    [("404", grzexc.MissingObjectError), ("403", grzexc.DownloadError)],
)
def test_head_object_sorts_the_head_error_if_the_get_request_succeeds(head_status: str, expected: type[Exception]):
    """The object may appear between the two requests. The GET response is then closed unread."""
    s3_client = _FailingClient(head_status)

    with pytest.raises(expected) as excinfo:
        head_object(s3_client, "bucket", "key")

    assert type(excinfo.value) is expected
    assert s3_client.body.closed


@pytest.mark.parametrize(
    ("head_status", "get_code", "expected"),
    [
        ("404", "NoSuchKey", grzexc.MissingSubmissionFileError),
        ("404", "NoSuchBucket", grzexc.ConfigurationError),
        ("403", "InvalidAccessKeyId", grzexc.ConfigurationError),
    ],
)
def test_get_metadata_upload_timestamp_sorts_the_error_of_the_get_request(
    head_status: str, get_code: str, expected: type[Exception]
):
    with pytest.raises(expected):
        get_metadata_upload_timestamp(_FailingClient(head_status, get_code), "bucket", "submission")
