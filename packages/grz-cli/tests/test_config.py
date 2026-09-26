"""Tests for ``KeyModel``, the keys section of the grz-cli config, which moved here from grz-common."""

import importlib.util
from pathlib import Path

import pytest
from grz_cli.models.config import EncryptConfig, KeyModel
from pydantic import ValidationError

PUBLIC_KEY = "-----BEGIN CRYPT4GH PUBLIC KEY-----\n7JZ9eRjhOo1zB8HfoQK1ULCR3Wpnl91hF2K8FtpmeQ8=\n-----END CRYPT4GH PUBLIC KEY-----\n"
"""A crypt4gh public key as text. ``Crypt4GHPublicKey`` only checks for the markers, so it need not decode."""

SUBMITTER_PRIVATE_KEY = "submitter-private-key"
"""Stands for a private key given inline. The model does not parse it, so it need not be a key."""


def test_grz_common_keeps_no_copy_of_the_key_models():
    assert importlib.util.find_spec("grz_common.models.keys") is None


def test_encrypt_config_reads_the_keys_section(tmp_path: Path):
    public_key_path = tmp_path / "grz.pub"
    public_key_path.write_text(PUBLIC_KEY)

    config = EncryptConfig.model_validate({"keys": {"grz_public_key_path": str(public_key_path)}})

    assert config.keys.grz_public_key_path == public_key_path


def test_a_config_with_the_removed_grz_private_key_path_still_loads():
    """grz-cli never read ``keys.grz_private_key_path``. Old configs still set it, so the model ignores it."""
    keys = KeyModel.model_validate({"grz_public_key": PUBLIC_KEY, "grz_private_key_path": "/no/such/key.sec"})

    assert "grz_private_key_path" not in keys.model_dump()


def test_neither_grz_public_key_nor_grz_public_key_path_fails():
    with pytest.raises(ValidationError, match="Either grz_public_key or grz_public_key_path must be set"):
        KeyModel()


def test_both_grz_public_key_and_grz_public_key_path_fails(tmp_path: Path):
    public_key_path = tmp_path / "grz.pub"
    public_key_path.write_text(PUBLIC_KEY)

    with pytest.raises(ValidationError, match="Only one of grz_public_key or grz_public_key_path must be set"):
        KeyModel(grz_public_key=PUBLIC_KEY, grz_public_key_path=public_key_path)


@pytest.mark.parametrize(
    "public_key",
    [
        PUBLIC_KEY,
        "-----BEGIN CRYPT4GH PUBLIC KEY-----\n7JZ9eRjhOo1zB8HfoQK1ULCR3Wpnl91hF2K8FtpmeQ8=\n",
        "7JZ9eRjhOo1zB8HfoQK1ULCR3Wpnl91hF2K8FtpmeQ8=\n-----END CRYPT4GH PUBLIC KEY-----\n",
    ],
    ids=["both markers", "BEGIN marker only", "END marker only"],
)
def test_grz_public_key_with_a_marker_passes(public_key: str):
    assert KeyModel(grz_public_key=public_key).grz_public_key == public_key


def test_grz_public_key_without_markers_fails():
    with pytest.raises(ValidationError, match="Invalid public key format"):
        KeyModel(grz_public_key="7JZ9eRjhOo1zB8HfoQK1ULCR3Wpnl91hF2K8FtpmeQ8=")


def test_both_submitter_private_key_and_submitter_private_key_path_fails(tmp_path: Path):
    private_key_path = tmp_path / "submitter.sec"
    private_key_path.write_text(SUBMITTER_PRIVATE_KEY)

    with pytest.raises(
        ValidationError, match="Only one of submitter_private_key or submitter_private_key_path must be set"
    ):
        KeyModel(
            grz_public_key=PUBLIC_KEY,
            submitter_private_key=SUBMITTER_PRIVATE_KEY,
            submitter_private_key_path=private_key_path,
        )


def test_the_submitter_private_key_is_masked_in_dumps():
    keys = KeyModel(grz_public_key=PUBLIC_KEY, submitter_private_key=SUBMITTER_PRIVATE_KEY)

    assert SUBMITTER_PRIVATE_KEY not in repr(keys)
    assert SUBMITTER_PRIVATE_KEY not in keys.model_dump_json()


def test_the_submitter_private_key_passphrase_is_masked_in_dumps():
    passphrase = "submitter-passphrase"
    keys = KeyModel(grz_public_key=PUBLIC_KEY, submitter_private_key_passphrase=passphrase)

    assert passphrase not in repr(keys)
    assert passphrase not in keys.model_dump_json()


def test_a_config_error_does_not_show_the_submitter_private_key(tmp_path: Path):
    private_key_path = tmp_path / "submitter.sec"
    private_key_path.write_text(SUBMITTER_PRIVATE_KEY)
    keys = {
        "grz_public_key": PUBLIC_KEY,
        "submitter_private_key": SUBMITTER_PRIVATE_KEY,
        "submitter_private_key_path": str(private_key_path),
    }

    with pytest.raises(ValidationError) as exc_info:
        EncryptConfig.model_validate({"keys": keys})

    assert SUBMITTER_PRIVATE_KEY not in str(exc_info.value)
