"""grzctl dump-config prints the loaded configuration as YAML, with secrets masked unless asked otherwise."""

import json
from pathlib import Path

import click.testing
import grzctl.cli
import pytest
import yaml
from grzctl.models.config import GrzctlConfig

PASSPHRASE = "author-passphrase"
CLIENT_SECRET = "pruefbericht-client-secret"
INBOX_S3_SECRET = "inbox-s3-secret"
ARCHIVE_S3_SECRET = "archive-s3-secret"


def _dump_config(config_path: Path, *args: str) -> str:
    runner = click.testing.CliRunner()
    result = runner.invoke(grzctl.cli.build_cli(), ["--config", str(config_path), "dump-config", *args])
    assert result.exit_code == 0, result.stderr
    return result.stdout


@pytest.fixture
def config_with_secrets_path(tmp_path: Path, offline_config: GrzctlConfig) -> Path:
    """A YAML config file with secrets in plain text: in an IgnoringBaseSettings, in an IgnoringBaseModel, and in S3."""
    data = offline_config.model_dump(mode="json", exclude_none=True)
    data["db"]["author"]["private_key_passphrase"] = PASSPHRASE
    data["pruefbericht"]["client_secret"] = CLIENT_SECRET
    data["leistungserbringer"]["000000000"]["inbox_buckets"]["inbox"]["secret"] = INBOX_S3_SECRET
    data["archives"]["consented"]["s3"]["secret"] = ARCHIVE_S3_SECRET

    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(data))
    return config_path


def test_dump_config_masks_secrets_by_default(config_with_secrets_path: Path):
    dumped = yaml.safe_load(_dump_config(config_with_secrets_path))

    assert dumped["db"]["author"]["private_key_passphrase"] == "**********"
    assert dumped["pruefbericht"]["client_secret"] == "**********"
    assert dumped["leistungserbringer"]["000000000"]["inbox_buckets"]["inbox"]["secret"] == "**********"
    assert dumped["archives"]["consented"]["s3"]["secret"] == "**********"


def test_dump_config_reveal_secrets_roundtrips(tmp_path: Path, config_with_secrets_path: Path):
    """The --reveal-secrets output loads back as a config file, and dumping that again prints the same YAML."""
    first = _dump_config(config_with_secrets_path, "--reveal-secrets")
    dumped = yaml.safe_load(first)
    assert dumped["db"]["author"]["private_key_passphrase"] == PASSPHRASE
    assert dumped["pruefbericht"]["client_secret"] == CLIENT_SECRET

    reloaded_path = tmp_path / "reloaded.yaml"
    reloaded_path.write_text(first)
    assert _dump_config(reloaded_path, "--reveal-secrets") == first
    assert GrzctlConfig.from_path(reloaded_path) == GrzctlConfig.from_path(config_with_secrets_path)


def test_json_dump_with_reveal_secrets_roundtrips(config_with_secrets_path: Path):
    """With reveal_secrets, a JSON dump loads back into an equal config, secrets included."""
    config = GrzctlConfig.from_path(config_with_secrets_path)

    dumped = config.model_dump_json(context={"reveal_secrets": True})

    assert GrzctlConfig.from_configuration(json.loads(dumped)) == config
