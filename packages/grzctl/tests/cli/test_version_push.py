"""Tests for `grzctl version-push`."""

from collections.abc import Iterator
from pathlib import Path
from typing import Any

import boto3
import grzctl.cli
import pytest
from click.testing import CliRunner
from grz_common.models.version import VersionFile
from moto import mock_aws

from .conftest import _database_config, _write_config

BUCKET = "inbox"
REGION = "us-east-1"
SUBMITTER_ID = "000000000"
INBOX_NAME = "inbox"


@pytest.fixture
def s3_client_mock() -> Iterator[Any]:
    """A moto-backed S3 client with a pre-created inbox bucket."""
    with mock_aws():
        client = boto3.client("s3", region_name=REGION)
        client.create_bucket(Bucket=BUCKET)
        yield client


@pytest.fixture
def config_path(tmp_path: Path) -> Path:
    """A GrzctlConfig with a single LE/inbox resolving to the mocked bucket."""
    config = _database_config(tmp_path, "sqlite:///:memory:")
    return _write_config(tmp_path, config)


def test_version_push_publishes_bundled_file(s3_client_mock, config_path):
    """Publishes the bundled version.json to the resolved inbox bucket."""
    runner = CliRunner()
    cli = grzctl.cli.build_cli()
    result = runner.invoke(
        cli,
        ["--config", str(config_path), "version-push", "-s", SUBMITTER_ID, "-b", INBOX_NAME],
    )

    assert result.exit_code == 0, result.output

    uploaded = s3_client_mock.get_object(Bucket=BUCKET, Key="version.json")
    content = uploaded["Body"].read().decode("utf-8")
    assert content == VersionFile.read_bundled_text()


def test_version_push_unknown_inbox_fails(s3_client_mock, config_path):
    """Fails clearly when the submitter/inbox pair isn't configured."""
    runner = CliRunner()
    cli = grzctl.cli.build_cli()
    result = runner.invoke(
        cli,
        ["--config", str(config_path), "version-push", "-s", SUBMITTER_ID, "-b", "does-not-exist"],
    )

    assert result.exit_code != 0
