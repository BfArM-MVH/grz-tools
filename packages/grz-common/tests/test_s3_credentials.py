"""The S3 client signs with the configured credentials, and without any it leaves the search to boto3."""

import os

import boto3
import pytest
from botocore.exceptions import PartialCredentialsError
from grz_common.models.s3 import S3Options
from grz_common.transfer import init_s3_client

S3 = {"endpoint_url": "https://s3.example.org", "bucket": "test-bucket"}


@pytest.fixture(autouse=True)
def aws_variables(monkeypatch):
    """Set the AWS variables, and keep boto3 from the files of this machine and from its cached default session."""
    monkeypatch.setattr(boto3, "DEFAULT_SESSION", None)
    monkeypatch.setenv("AWS_SHARED_CREDENTIALS_FILE", os.devnull)
    monkeypatch.setenv("AWS_CONFIG_FILE", os.devnull)
    monkeypatch.delenv("AWS_PROFILE", raising=False)
    monkeypatch.delenv("AWS_SESSION_TOKEN", raising=False)
    monkeypatch.setenv("AWS_ACCESS_KEY_ID", "aws-key")
    monkeypatch.setenv("AWS_SECRET_ACCESS_KEY", "aws-secret")


def _signed_url(s3_options: S3Options) -> str:
    """Sign a URL offline. It names the access key, and signing fails without credentials.

    :param s3_options: Options of the S3 client.
    :returns: The signed URL.
    """
    return init_s3_client(s3_options).generate_presigned_url(
        "get_object", Params={"Bucket": s3_options.bucket, "Key": "key"}
    )


def test_without_configured_credentials_boto3_finds_the_aws_variables():
    assert "aws-key" in _signed_url(S3Options.model_validate(S3))


def test_configured_credentials_win_over_the_aws_variables():
    url = _signed_url(S3Options.model_validate({**S3, "access_key": "file-key", "secret": "file-secret"}))

    assert "file-key" in url
    assert "aws-key" not in url


@pytest.mark.parametrize("credential", [{"access_key": "file-key"}, {"secret": "file-secret"}], ids=["key", "secret"])
def test_one_configured_credential_does_not_fall_back_to_the_aws_variables(credential: dict[str, str]):
    with pytest.raises(PartialCredentialsError):
        init_s3_client(S3Options.model_validate({**S3, **credential}))
