"""Tests for the grzctl ``encrypt`` command."""

from base64 import b64decode
from pathlib import Path
from unittest.mock import MagicMock, patch

import click.testing
import crypt4gh.keys
import crypt4gh.keys.c4gh
import grz_common.exceptions as grzexc
import grzctl.cli
import pytest
import yaml

from .conftest import PRUEFBERICHT

SIGNING_KEY_PASSPHRASE = "signing-key-passphrase"


@pytest.fixture
def signing_key_path(tmp_path) -> Path:
    """A crypt4gh private key, encrypted with ``SIGNING_KEY_PASSPHRASE``."""
    private_key_path = tmp_path / "grz.sec"
    crypt4gh.keys.c4gh.generate(private_key_path, tmp_path / "grz.pub", SIGNING_KEY_PASSPHRASE.encode(), comment=None)
    return private_key_path


def _config(public_key: str, unread_file: str, signing_key: dict[str, str]) -> dict:
    return {
        "leistungserbringer": {"000000000": {"inbox_buckets": {"inbox": {"private_key_path": unread_file}}}},
        "archives": {
            "consented": {"s3": {"bucket": "consented"}, "public_key": public_key},
            "non_consented": {"s3": {"bucket": "non_consented"}, "public_key_path": unread_file},
            **signing_key,
        },
        "db": {"database_url": "sqlite:///:memory:", "author": {"name": "test"}},
        "pruefbericht": PRUEFBERICHT,
        "identifiers": {"grz": "GRZT00000"},
    }


def _write_config(tmp_path: Path, config: dict) -> Path:
    config_path = tmp_path / "config.yaml"
    with open(config_path, "w") as f:
        yaml.dump(config, f)
    return config_path


@pytest.fixture
def grzctl_config_path(tmp_path, signing_key_path, crypt4gh_public_key, unread_file):
    signing_key = {"signing_key_path": str(signing_key_path), "signing_key_passphrase": SIGNING_KEY_PASSPHRASE}
    return _write_config(tmp_path, _config(crypt4gh_public_key, unread_file, signing_key))


def _submission_dir(tmp_path: Path) -> Path:
    submission_dir = tmp_path / "submission"
    for sub in ("metadata", "files", "logs", "encrypted_files"):
        (submission_dir / sub).mkdir(parents=True)
    return submission_dir


def _invoke_encrypt(config_path: Path, submission_dir: Path, mock_worker_cls: MagicMock) -> click.testing.Result:
    """Run ``grzctl encrypt`` against a mocked ``Worker`` for a consented submission."""
    mock_worker = mock_worker_cls.return_value
    mock_submission = mock_worker.parse_submission.return_value
    mock_submission.metadata.content.submission_id = "S1"
    mock_submission.metadata.content.consents_to_research.return_value = True

    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()
    return runner.invoke(
        cli, ["--config", str(config_path), "encrypt", "--submission-dir", str(submission_dir), "--no-update-db"]
    )


def test_encrypt_uses_an_inline_archive_public_key(tmp_path, grzctl_config_path, crypt4gh_public_key):
    """The consented archive's public key may be given inline instead of as a file path."""
    with patch("grzctl.commands.encrypt.Worker") as mock_worker_cls:
        result = _invoke_encrypt(grzctl_config_path, _submission_dir(tmp_path), mock_worker_cls)

    assert result.exit_code == 0, result.output
    mock_worker_cls.return_value.encrypt.assert_called_once()
    recipient_public_key = mock_worker_cls.return_value.encrypt.call_args.kwargs["recipient_public_key"]
    assert recipient_public_key.public_bytes_raw() == b64decode(crypt4gh_public_key.splitlines()[1])


def test_encrypt_signs_with_the_signing_key(tmp_path, grzctl_config_path, signing_key_path, monkeypatch):
    """The signing key reaches ``Worker.encrypt`` loaded, decrypted with the configured passphrase."""
    monkeypatch.setenv("C4GH_PASSPHRASE", "wrong-passphrase")

    with patch("grzctl.commands.encrypt.Worker") as mock_worker_cls:
        result = _invoke_encrypt(grzctl_config_path, _submission_dir(tmp_path), mock_worker_cls)

    assert result.exit_code == 0, result.output
    encrypt_kwargs = mock_worker_cls.return_value.encrypt.call_args.kwargs
    expected_key = crypt4gh.keys.get_private_key(signing_key_path, lambda: SIGNING_KEY_PASSPHRASE)
    assert encrypt_kwargs["submitter_private_key"].private_bytes_raw() == expected_key


def test_encrypt_signs_with_an_inline_signing_key(tmp_path, signing_key_path, crypt4gh_public_key, unread_file):
    signing_key = {"signing_key": signing_key_path.read_text(), "signing_key_passphrase": SIGNING_KEY_PASSPHRASE}
    config_path = _write_config(tmp_path, _config(crypt4gh_public_key, unread_file, signing_key))

    with patch("grzctl.commands.encrypt.Worker") as mock_worker_cls:
        result = _invoke_encrypt(config_path, _submission_dir(tmp_path), mock_worker_cls)

    assert result.exit_code == 0, result.output
    expected_key = crypt4gh.keys.get_private_key(signing_key_path, lambda: SIGNING_KEY_PASSPHRASE)
    signing_key = mock_worker_cls.return_value.encrypt.call_args.kwargs["submitter_private_key"]
    assert signing_key.private_bytes_raw() == expected_key


def test_encrypt_fails_if_the_signing_key_cannot_be_loaded(tmp_path, crypt4gh_public_key, unread_file):
    """A signing key that cannot be loaded is a configuration error, which the GRZ has to fix."""
    not_a_key_path = tmp_path / "not_a_key.sec"
    not_a_key_path.write_text("not a key")
    config_path = _write_config(
        tmp_path, _config(crypt4gh_public_key, unread_file, {"signing_key_path": str(not_a_key_path)})
    )

    with patch("grzctl.commands.encrypt.Worker") as mock_worker_cls:
        result = _invoke_encrypt(config_path, _submission_dir(tmp_path), mock_worker_cls)

    assert isinstance(result.exception, grzexc.ConfigurationError), result.output
    assert f"Secret key {not_a_key_path} cannot be read" in str(result.exception)
    mock_worker_cls.return_value.encrypt.assert_not_called()
