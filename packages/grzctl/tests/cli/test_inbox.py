"""Tests for `grzctl inbox push-version`."""

from collections.abc import Iterator
from pathlib import Path
from typing import Any

import boto3
import cryptography.hazmat.primitives.serialization as cryptser
import grzctl.cli
import pytest
from click.testing import CliRunner
from cryptography.hazmat.primitives.asymmetric.ed25519 import Ed25519PrivateKey
from grz_common.models.version import VERSION_FILE_KEY, VersionFile
from grzctl.models.config import GrzctlConfig
from moto import mock_aws

from .conftest import _grzctl_archives, _write_config

REGION = "us-east-1"
BUCKET_A = "inbox-a"
BUCKET_B = "inbox-b"


def _two_inbox_config(tmp_path: Path) -> GrzctlConfig:
    """A config with two LEs, each with one inbox, resolving to distinct buckets."""
    private_key = Ed25519PrivateKey.generate()
    private_key_path = tmp_path / "grz.sec"
    with open(private_key_path, "wb") as private_key_file:
        private_key_file.write(
            private_key.private_bytes(
                encoding=cryptser.Encoding.PEM,
                format=cryptser.PrivateFormat.OpenSSH,
                encryption_algorithm=cryptser.NoEncryption(),
            )
        )

    public_key = private_key.public_key()
    public_key_path = tmp_path / "grz.pub"
    with open(public_key_path, "wb") as public_key_file:
        public_key_file.write(
            public_key.public_bytes(encoding=cryptser.Encoding.OpenSSH, format=cryptser.PublicFormat.OpenSSH)
        )

    inbox = {"private_key_path": str(private_key_path.resolve())}
    return GrzctlConfig(
        leistungserbringer={
            "000000000": {"inbox_buckets": {"inbox": {"bucket": BUCKET_A, **inbox}}},
            "111111111": {"inbox_buckets": {"inbox": {"bucket": BUCKET_B, **inbox}}},
        },
        archives=_grzctl_archives(
            public_key_path=str(public_key_path.resolve()),
            signing_key_path=str(private_key_path.resolve()),
        ),
        db={
            "database_url": "sqlite:///:memory:",
            "author": {
                "name": "alice",
                "private_key_path": str(private_key_path.resolve()),
                "private_key_passphrase": "",
            },
            "known_public_keys_file": str(public_key_path.resolve()),
        },
        pruefbericht={},
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
    return _write_config(tmp_path, _two_inbox_config(tmp_path))


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
