"""Tests for `grzctl inbox push-version`."""

from collections.abc import Iterator
from pathlib import Path
from typing import Any

import boto3
import grzctl.cli
import pytest
from click.testing import CliRunner
from grz_common.models.version import VERSION_FILE_KEY, VersionFile
from grzctl.models.config import GrzctlConfig
from moto import mock_aws

from .conftest import _grzctl_archives, _write_config

REGION = "us-east-1"
BUCKET_A = "inbox-a"
BUCKET_B = "inbox-b"


def _two_inbox_config() -> GrzctlConfig:
    """A config with two LEs, each with one inbox, resolving to distinct buckets."""
    return GrzctlConfig(
        leistungserbringer={
            "000000000": {"inbox_buckets": {"inbox": {"bucket": BUCKET_A, "private_key_path": "unused"}}},
            "111111111": {"inbox_buckets": {"inbox": {"bucket": BUCKET_B, "private_key_path": "unused"}}},
        },
        archives=_grzctl_archives(),
        db={"database_url": "sqlite:///:memory:", "author": {"name": "alice"}},
        pruefbericht={},
        keys={"grz_private_key_path": "unused"},
        identifiers={"grz": "GRZK00007"},
    )


@pytest.fixture
def s3_client_mock() -> Iterator[Any]:
    """A moto-backed S3 client with both inbox buckets pre-created."""
    with mock_aws():
        client = boto3.client("s3", region_name=REGION)
        client.create_bucket(Bucket=BUCKET_A)
        client.create_bucket(Bucket=BUCKET_B)
        yield client


@pytest.fixture
def config_path(tmp_path: Path) -> Path:
    return _write_config(tmp_path, _two_inbox_config())


def test_push_version_publishes_to_every_inbox(s3_client_mock, config_path):
    """Publishes the bundled version.json to every configured LE/inbox, not just one."""
    runner = CliRunner()
    cli = grzctl.cli.build_cli()
    result = runner.invoke(cli, ["--config", str(config_path), "inbox", "push-version"])

    assert result.exit_code == 0, result.output

    bundled = VersionFile.read_bundled_text()
    for bucket in (BUCKET_A, BUCKET_B):
        uploaded = s3_client_mock.get_object(Bucket=bucket, Key=VERSION_FILE_KEY)
        assert uploaded["Body"].read().decode("utf-8") == bundled


def test_push_version_reports_failure_for_missing_bucket(config_path):
    """A missing/unreachable bucket is reported and fails the command, without a real S3 backend."""
    with mock_aws():
        runner = CliRunner()
        cli = grzctl.cli.build_cli()
        result = runner.invoke(cli, ["--config", str(config_path), "inbox", "push-version"])

    assert result.exit_code != 0
