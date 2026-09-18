"""An error of the S3 client becomes the failure it stands for."""

import pytest
from botocore.exceptions import ClientError, EndpointConnectionError, NoCredentialsError
from grz_common.exceptions import ConfigurationError, DownloadError, MissingObjectError, UploadError
from grz_common.pipeline.components.s3 import head_object, s3_error


def _client_error(code: str) -> ClientError:
    return ClientError({"Error": {"Code": code, "Message": code}}, "HeadObject")


class _FailingClient:
    def __init__(self, error: Exception):
        self.error = error

    def head_object(self, **_kwargs):
        raise self.error


@pytest.mark.parametrize("code", ["InvalidAccessKeyId", "SignatureDoesNotMatch", "NoSuchBucket"])
def test_an_error_that_only_a_faulty_setup_causes_is_a_configuration_error(code: str):
    assert isinstance(s3_error(_client_error(code), "Reading s3://bucket/key"), ConfigurationError)


def test_missing_credentials_are_a_configuration_error():
    assert isinstance(s3_error(NoCredentialsError(), "Reading s3://bucket/key"), ConfigurationError)


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
    failure = s3_error(error, "Upload to s3://bucket/key", UploadError)

    assert type(failure) is UploadError
    assert str(failure).startswith("Upload to s3://bucket/key failed: ")


def test_head_object_reports_a_missing_object():
    with pytest.raises(MissingObjectError):
        head_object(_FailingClient(_client_error("404")), "bucket", "key")


def test_head_object_reports_a_refused_read_as_a_failed_download():
    with pytest.raises(DownloadError) as excinfo:
        head_object(_FailingClient(_client_error("AccessDenied")), "bucket", "key")

    assert not isinstance(excinfo.value, MissingObjectError)
