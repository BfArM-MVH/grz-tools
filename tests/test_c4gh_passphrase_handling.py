"""Unit tests for Crypt4GH utility wrapper."""

import os
from stat import S_IRUSR, S_IWUSR
from unittest.mock import patch

import crypt4gh.keys.c4gh
import pytest
from grz_common.utils.crypt import Crypt4GH


@pytest.fixture
def encrypted_dummy_key(tmp_path) -> tuple[str, bytes]:
    """Generates a Crypt4GH keypair encrypted with a passphrase."""
    sec_key_path = tmp_path / "dummy_encrypted.sec"
    pub_key_path = tmp_path / "dummy_encrypted.pub"
    passphrase = b"my-secret-test-passphrase"

    crypt4gh.keys.c4gh.generate(str(sec_key_path), str(pub_key_path), passphrase=passphrase)
    os.chmod(str(sec_key_path), S_IRUSR | S_IWUSR)
    return str(sec_key_path), passphrase


def test_retrieve_private_key_with_explicit_passphrase(encrypted_dummy_key, monkeypatch):
    sec_key_path, passphrase = encrypted_dummy_key
    monkeypatch.delenv("C4GH_PASSPHRASE", raising=False)

    key_bytes = Crypt4GH.retrieve_private_key(sec_key_path, passphrase=passphrase.decode("utf-8"))

    assert isinstance(key_bytes, bytes)
    assert len(key_bytes) == 32


def test_retrieve_private_key_with_envvar_fallback(encrypted_dummy_key, monkeypatch):
    sec_key_path, passphrase = encrypted_dummy_key
    monkeypatch.setenv("C4GH_PASSPHRASE", passphrase.decode("utf-8"))

    key_bytes = Crypt4GH.retrieve_private_key(sec_key_path, passphrase=None)

    assert isinstance(key_bytes, bytes)
    assert len(key_bytes) == 32


def test_retrieve_private_key_missing_passphrase(encrypted_dummy_key, monkeypatch):
    sec_key_path, _ = encrypted_dummy_key
    monkeypatch.delenv("C4GH_PASSPHRASE", raising=False)

    with patch("grz_common.utils.crypt.getpass") as mock_getpass:
        mock_getpass.side_effect = RuntimeError("Interactive prompt triggered in CI")

        with pytest.raises(RuntimeError, match="Interactive prompt triggered in CI"):
            Crypt4GH.retrieve_private_key(sec_key_path, passphrase=None)

        mock_getpass.assert_called_once()
