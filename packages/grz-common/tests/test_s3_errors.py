"""An error of the S3 client becomes the failure it stands for."""

import datetime
import io

import boto3
import grz_common.exceptions as grzexc
import pytest
from botocore.exceptions import ClientError, EndpointConnectionError, NoCredentialsError
from grz_common.transfer import get_metadata_upload_timestamp, head_object, s3_errors


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


def _failure_of(error: Exception, transfer_error: type[grzexc.TransferError] = grzexc.TransferError) -> Exception:
    """Return what ``s3_errors`` raises for ``error``."""
    with pytest.raises(grzexc.GrzError) as excinfo, s3_errors("Upload to s3://bucket/key", transfer_error):
        raise error
    return excinfo.value


@pytest.mark.parametrize("code", ["InvalidAccessKeyId", "SignatureDoesNotMatch", "NoSuchBucket"])
def test_an_error_that_only_a_faulty_setup_causes_is_a_configuration_error(code: str):
    assert isinstance(_failure_of(_client_error(code)), grzexc.ConfigurationError)


def test_missing_credentials_are_a_configuration_error():
    assert isinstance(_failure_of(NoCredentialsError()), grzexc.ConfigurationError)


def test_a_bucket_name_that_botocore_rejects_is_a_configuration_error():
    """The S3 client checks the bucket name before it sends a request, so no S3 server is needed."""
    s3_client = boto3.client("s3", region_name="us-east-1", aws_access_key_id="key", aws_secret_access_key="secret")

    with pytest.raises(grzexc.ConfigurationError):
        head_object(s3_client, "not a bucket name", "key")


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
    failure = _failure_of(error, grzexc.UploadError)

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


@pytest.mark.parametrize("head_status", ["404", "403"])
def test_head_object_reports_a_failed_read_if_the_get_request_succeeds(head_status: str):
    """The object may appear between the two requests. The GET response is then closed unread."""
    s3_client = _FailingClient(head_status)

    with pytest.raises(grzexc.DownloadError) as excinfo:
        head_object(s3_client, "bucket", "key")

    assert not isinstance(excinfo.value, grzexc.MissingObjectError)
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


UPLOADED = datetime.datetime(2024, 7, 15, tzinfo=datetime.UTC)


class _InboxClient:
    """Answer a HEAD request for metadata.json with ``content_length``, and list the keys in ``keys``."""

    def __init__(self, content_length: int = 2, keys: tuple[str, ...] = ()):
        self.content_length = content_length
        self.keys = keys

    def head_object(self, **_kwargs):
        return {"ContentLength": self.content_length, "LastModified": UPLOADED}

    def list_objects_v2(self, Bucket: str, Prefix: str):  # noqa: N803
        contents = [{"Key": key} for key in self.keys if key.startswith(Prefix)]
        # S3 leaves out Contents when no key matches
        return {"Contents": contents} if contents else {}


def test_get_metadata_upload_timestamp_returns_the_time_of_upload():
    assert get_metadata_upload_timestamp(_InboxClient(keys=("submission/version",)), "bucket", "submission") == UPLOADED


@pytest.mark.parametrize(
    "inbox",
    [
        _InboxClient(keys=("submission/cleaning",)),
        _InboxClient(content_length=0, keys=("submission/cleaned",)),
        _InboxClient(content_length=0),
    ],
    ids=["being-cleaned", "cleaned", "emptied"],
)
def test_get_metadata_upload_timestamp_refuses_a_cleaned_submission(inbox: _InboxClient):
    """An emptied metadata.json carries the time of cleaning, not the time of upload."""
    with pytest.raises(grzexc.SubmissionCleanedError):
        get_metadata_upload_timestamp(inbox, "bucket", "submission")
