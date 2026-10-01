"""Unit tests for Crypt4GH utility wrapper."""

import os
from stat import S_IRUSR, S_IWUSR

import crypt4gh.keys.c4gh
import grz_common.exceptions as grzexc
import pytest
from cryptography.hazmat.primitives.asymmetric.x25519 import X25519PrivateKey
from grz_common.utils.crypt import Crypt4GH


@pytest.fixture
def encrypted_dummy_key(tmp_path) -> tuple[str, bytes]:
    """Generates a Crypt4GH keypair encrypted with a passphrase."""
    sec_key_path = tmp_path / "dummy_encrypted.sec"
    pub_key_path = tmp_path / "dummy_encrypted.pub"
    passphrase = b"my-secret-test-passphrase"

    # c4gh modifies umask _process wide_ so we have to be able to undo that…
    prev_umask = os.umask(0)
    os.umask(prev_umask)

    try:
        crypt4gh.keys.c4gh.generate(str(sec_key_path), str(pub_key_path), passphrase=passphrase, comment=b"dummy key")
    finally:
        os.umask(prev_umask)

    os.chmod(str(sec_key_path), S_IRUSR | S_IWUSR)
    return str(sec_key_path), passphrase


def test_retrieve_private_key_with_envvar(encrypted_dummy_key, monkeypatch):
    sec_key_path, passphrase = encrypted_dummy_key
    monkeypatch.setenv("C4GH_PASSPHRASE", passphrase.decode("utf-8"))

    private_key = Crypt4GH.retrieve_private_key(sec_key_path)

    assert isinstance(private_key, X25519PrivateKey)


def test_a_missing_private_key_is_a_configuration_error(tmp_path):
    with pytest.raises(grzexc.ConfigurationError):
        Crypt4GH.retrieve_private_key(tmp_path / "missing.sec")


def test_a_wrong_passphrase_is_a_configuration_error(encrypted_dummy_key, monkeypatch):
    """crypt4gh exits the process for a wrong passphrase, which would pass every ``except Exception``."""
    sec_key_path, _ = encrypted_dummy_key
    monkeypatch.setenv("C4GH_PASSPHRASE", "not the passphrase")

    with pytest.raises(grzexc.ConfigurationError):
        Crypt4GH.retrieve_private_key(sec_key_path)


def test_a_missing_public_key_is_a_configuration_error(tmp_path):
    with pytest.raises(grzexc.ConfigurationError):
        Crypt4GH.retrieve_public_key(tmp_path / "missing.pub")
