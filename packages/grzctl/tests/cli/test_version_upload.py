"""Tests for `grzctl version-upload`."""

import json
from collections.abc import Iterator
from datetime import UTC, datetime
from typing import Any

import boto3
import grzctl.cli
import pytest
import responses
from click.testing import CliRunner
from grz_common.models.s3 import S3ConfigModel
from grzctl.commands.version_upload import VERSION_FILE_URL
from moto import mock_aws

BUCKET = "test-inbox-bucket"
REGION = "us-east-1"

VERSION_CONTENT = {
    "schema_version": 1,
    "grzcli_version": [
        {
            "minimal_version": "1.5.0",
            "recommended_version": "1.7.0",
            "max_version": "2.0.0",
            "enforced_from": datetime.now(UTC).isoformat(),
        }
    ],
    "metadata_version": [],
}


@pytest.fixture
def requests_mock(assert_all_requests_are_fired: bool = False):
    with responses.RequestsMock(assert_all_requests_are_fired=assert_all_requests_are_fired) as rsps:
        yield rsps


@pytest.fixture
def s3_client_mock() -> Iterator[Any]:
    """A moto-backed S3 client with a pre-created inbox bucket."""
    with mock_aws():
        client = boto3.client("s3", region_name=REGION)
        client.create_bucket(Bucket=BUCKET)
        yield client


@pytest.fixture
def config_file(tmp_path):
    config = S3ConfigModel(
        s3={
            "endpoint_url": "https://s3.amazonaws.com",
            "bucket": BUCKET,
            "access_key": "testing",
            "secret": "testing",
        }
    )
    path = tmp_path / "config-inbox.yaml"
    with open(path, "w", encoding="utf-8") as fd:
        config.to_yaml(fd)
    return path


def test_version_upload_publishes_downloaded_file(s3_client_mock, config_file, requests_mock):
    """Downloads version.json from GitHub and uploads it to the configured bucket."""
    requests_mock.get(VERSION_FILE_URL, json=VERSION_CONTENT)

    runner = CliRunner()
    cli = grzctl.cli.build_cli()
    result = runner.invoke(
        cli,
        ["version-upload", "--config-file", str(config_file)],
    )

    assert result.exit_code == 0, result.output

    uploaded = s3_client_mock.get_object(Bucket=BUCKET, Key="version.json")
    assert json.loads(uploaded["Body"].read().decode("utf-8")) == VERSION_CONTENT


def test_version_upload_rejects_invalid_policy(s3_client_mock, config_file, requests_mock):
    """An invalid version.json from GitHub is not uploaded."""
    requests_mock.get(
        VERSION_FILE_URL,
        json={"schema_version": 1, "grzcli_version": [{"minimal_version": "not-a-version"}]},
    )

    runner = CliRunner()
    cli = grzctl.cli.build_cli()
    result = runner.invoke(
        cli,
        ["version-upload", "--config-file", str(config_file)],
    )

    assert result.exit_code != 0

    with pytest.raises(s3_client_mock.exceptions.NoSuchKey):
        s3_client_mock.get_object(Bucket=BUCKET, Key="version.json")
