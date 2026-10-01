"""Tests for loading crypt4gh keys in memory, the passphrase order of private keys, and encrypting and
decrypting with the loaded keys.
"""

from pathlib import Path
from unittest.mock import patch

import crypt4gh.header
import crypt4gh.keys
import crypt4gh.keys.c4gh
import cryptography.hazmat.primitives.serialization as cryptser
import grz_common.exceptions as grzexc
import pytest
from cryptography.hazmat.primitives.asymmetric.ed25519 import Ed25519PrivateKey
from cryptography.hazmat.primitives.asymmetric.x25519 import X25519PrivateKey
from grz_common.utils.crypt import Crypt4GH

PASSPHRASE = "right-passphrase"


def _generate_key_pair(tmp_path: Path, name: str, passphrase: str | None = None) -> tuple[Path, Path]:
    """Write a crypt4gh key pair, with the private key encrypted by *passphrase* if given."""
    private_key_path = tmp_path / f"{name}.sec"
    public_key_path = tmp_path / f"{name}.pub"
    crypt4gh.keys.c4gh.generate(
        private_key_path, public_key_path, passphrase.encode() if passphrase else None, comment=None
    )
    return private_key_path, public_key_path


@pytest.fixture
def plain_key_pair(tmp_path: Path) -> tuple[Path, Path]:
    return _generate_key_pair(tmp_path, "plain")


@pytest.fixture
def encrypted_key_pair(tmp_path: Path) -> tuple[Path, Path]:
    return _generate_key_pair(tmp_path, "encrypted", PASSPHRASE)


@pytest.fixture
def no_prompt(monkeypatch):
    """Fail the test if the passphrase prompt opens, and unset ``C4GH_PASSPHRASE``."""
    monkeypatch.delenv("C4GH_PASSPHRASE", raising=False)

    def _fail(*args, **kwargs):
        raise AssertionError("the passphrase prompt must not open")

    monkeypatch.setattr("grz_common.utils.crypt.getpass", _fail)


def test_load_private_key_matches_crypt4gh_for_a_crypt4gh_key(plain_key_pair, no_prompt):
    private_key_path, _ = plain_key_pair

    loaded = Crypt4GH.load_private_key(private_key_path.read_text())

    assert loaded.private_bytes_raw() == crypt4gh.keys.get_private_key(private_key_path, None)


@pytest.mark.parametrize("passphrase", [None, PASSPHRASE], ids=["unencrypted", "encrypted"])
def test_load_private_key_matches_crypt4gh_for_an_openssh_key(tmp_path: Path, no_prompt, passphrase: str | None):
    private_key_path = tmp_path / "id_ed25519"
    private_key_path.write_bytes(
        Ed25519PrivateKey.generate().private_bytes(
            encoding=cryptser.Encoding.PEM,
            format=cryptser.PrivateFormat.OpenSSH,
            encryption_algorithm=cryptser.BestAvailableEncryption(passphrase.encode())
            if passphrase
            else cryptser.NoEncryption(),
        )
    )

    loaded = Crypt4GH.load_private_key(private_key_path.read_text(), passphrase=passphrase)

    assert loaded.private_bytes_raw() == crypt4gh.keys.get_private_key(private_key_path, lambda: passphrase)


def test_load_private_key_writes_no_file(plain_key_pair, no_prompt):
    """The key is parsed from its text, so it never lands on disk."""
    private_key_text = plain_key_pair[0].read_text()

    with (
        patch("builtins.open", side_effect=AssertionError("no file may be opened")),
        patch("tempfile.NamedTemporaryFile", side_effect=AssertionError("no temporary file may be written")),
        patch("tempfile.mkstemp", side_effect=AssertionError("no temporary file may be written")),
    ):
        loaded = Crypt4GH.load_private_key(private_key_text)

    assert isinstance(loaded, X25519PrivateKey)


def test_configured_passphrase_comes_before_the_environment(encrypted_key_pair, monkeypatch, no_prompt):
    monkeypatch.setenv("C4GH_PASSPHRASE", "wrong-passphrase")
    private_key_path, _ = encrypted_key_pair

    loaded = Crypt4GH.load_private_key(private_key_path.read_text(), passphrase=PASSPHRASE)

    assert loaded.private_bytes_raw() == crypt4gh.keys.get_private_key(private_key_path, lambda: PASSPHRASE)


def test_environment_passphrase_comes_before_the_prompt(encrypted_key_pair, no_prompt, monkeypatch):
    monkeypatch.setenv("C4GH_PASSPHRASE", PASSPHRASE)
    private_key_path, _ = encrypted_key_pair

    loaded = Crypt4GH.load_private_key(private_key_path.read_text())

    assert loaded.private_bytes_raw() == crypt4gh.keys.get_private_key(private_key_path, lambda: PASSPHRASE)


def test_prompt_asks_for_the_passphrase_last(encrypted_key_pair, monkeypatch):
    monkeypatch.delenv("C4GH_PASSPHRASE", raising=False)
    prompts = []

    def _getpass(prompt):
        prompts.append(prompt)
        return PASSPHRASE

    monkeypatch.setattr("grz_common.utils.crypt.getpass", _getpass)
    private_key_path, _ = encrypted_key_pair

    Crypt4GH.load_private_key(private_key_path.read_text(), key_name="inbox.private_key")

    assert prompts == ["Passphrase for inbox.private_key: "]


def test_retrieve_private_key_uses_the_configured_passphrase(encrypted_key_pair, monkeypatch, no_prompt):
    monkeypatch.setenv("C4GH_PASSPHRASE", "wrong-passphrase")
    private_key_path, _ = encrypted_key_pair

    loaded = Crypt4GH.retrieve_private_key(private_key_path, passphrase=PASSPHRASE)

    assert loaded.private_bytes_raw() == crypt4gh.keys.get_private_key(private_key_path, lambda: PASSPHRASE)


def test_wrong_passphrase_raises_instead_of_exiting(encrypted_key_pair, no_prompt):
    """crypt4gh exits the process on a wrong passphrase, which would skip the caller's error handling."""
    private_key_path, _ = encrypted_key_pair

    with pytest.raises(grzexc.ConfigurationError, match=r"Secret key inbox\.private_key cannot be read with the given"):
        Crypt4GH.load_private_key(private_key_path.read_text(), passphrase="wrong", key_name="inbox.private_key")


@pytest.mark.parametrize(
    "private_key",
    ["not a key", "-----BEGIN SOMETHING-----\nbm90IGEga2V5\n-----END SOMETHING-----\n"],
)
def test_unsupported_key_fails(private_key, no_prompt):
    with pytest.raises(grzexc.ConfigurationError, match=r"Secret key \(inline\) cannot be read"):
        Crypt4GH.load_private_key(private_key)


def test_load_public_key_matches_crypt4gh_for_a_crypt4gh_key(plain_key_pair):
    _, public_key_path = plain_key_pair

    loaded = Crypt4GH.load_public_key(public_key_path.read_text())

    assert loaded.public_bytes_raw() == crypt4gh.keys.get_public_key(public_key_path)


def test_load_public_key_matches_crypt4gh_for_an_openssh_key(tmp_path: Path):
    public_key_path = tmp_path / "id_ed25519.pub"
    public_key_path.write_bytes(
        Ed25519PrivateKey.generate()
        .public_key()
        .public_bytes(encoding=cryptser.Encoding.OpenSSH, format=cryptser.PublicFormat.OpenSSH)
        + b" comment\n"
    )

    loaded = Crypt4GH.load_public_key(public_key_path.read_text())

    assert loaded.public_bytes_raw() == crypt4gh.keys.get_public_key(public_key_path)


def test_retrieve_public_key_matches_crypt4gh(plain_key_pair):
    _, public_key_path = plain_key_pair

    loaded = Crypt4GH.retrieve_public_key(public_key_path)

    assert loaded.public_bytes_raw() == crypt4gh.keys.get_public_key(public_key_path)


@pytest.mark.parametrize(
    "public_key",
    [
        "not a key",
        "-----BEGIN CRYPT4GH PUBLIC KEY-----\nbm90IGEga2V5\n-----END CRYPT4GH PUBLIC KEY-----\n",
        "ssh-ed25519 bm90IGEga2V5",
        "ssh-rsa bm90IGEga2V5",
    ],
    ids=["no key", "not 32 bytes", "malformed ssh-ed25519", "ssh-rsa"],
)
def test_unsupported_public_key_fails(public_key):
    with pytest.raises(grzexc.ConfigurationError, match=r"Public key \(inline\) cannot be read"):
        Crypt4GH.load_public_key(public_key)


def test_a_private_key_is_no_public_key(plain_key_pair):
    """The error names the key, but does not show it."""
    private_key_text = plain_key_pair[0].read_text()

    with pytest.raises(grzexc.ConfigurationError, match=r"Public key keys\.grz_public_key cannot be read") as exc_info:
        Crypt4GH.load_public_key(private_key_text, key_name="keys.grz_public_key")

    for line in private_key_text.splitlines():
        assert line not in str(exc_info.value)


def _encrypt(tmp_path: Path, public_key_path: Path, sender_private_key: X25519PrivateKey | None = None) -> Path:
    encrypted_path = tmp_path / "file.txt.c4gh"
    plain_path = tmp_path / "file.txt"
    plain_path.write_text("content")
    keys = Crypt4GH.prepare_c4gh_keys(Crypt4GH.retrieve_public_key(public_key_path), sender_private_key)
    Crypt4GH.encrypt_file(plain_path, encrypted_path, keys)
    return encrypted_path


def test_decrypt_file_with_a_key_that_does_not_open_the_header_fails(tmp_path: Path, plain_key_pair, no_prompt):
    """A key that does not fit is named as the reason, since crypt4gh's own message does not say so."""
    other_private_key_path, _ = _generate_key_pair(tmp_path, "other")
    encrypted_path = _encrypt(tmp_path, plain_key_pair[1])

    with pytest.raises(grzexc.DecryptionError, match="the private key does not open its Crypt4GH header"):
        Crypt4GH.decrypt_file(
            encrypted_path, tmp_path / "decrypted.txt", Crypt4GH.retrieve_private_key(other_private_key_path)
        )


def test_decrypt_file_fails_for_a_file_that_is_not_crypt4gh(tmp_path: Path, plain_key_pair, no_prompt):
    not_encrypted_path = tmp_path / "file.txt"
    not_encrypted_path.write_text("not a crypt4gh file")

    with pytest.raises(grzexc.DecryptionError, match="Not a CRYPT4GH formatted file"):
        Crypt4GH.decrypt_file(
            not_encrypted_path, tmp_path / "decrypted.txt", Crypt4GH.retrieve_private_key(plain_key_pair[0])
        )


def test_prepare_c4gh_keys_signs_with_the_sender_private_key(tmp_path: Path, plain_key_pair, no_prompt):
    """The header names the sender's public key, so the recipient can check who encrypted the file."""
    _, public_key_path = plain_key_pair
    sender_private_key_path, sender_public_key_path = _generate_key_pair(tmp_path, "sender")

    encrypted_path = _encrypt(tmp_path, public_key_path, Crypt4GH.retrieve_private_key(sender_private_key_path))

    with open(encrypted_path, "rb") as encrypted_file:
        (packet,) = crypt4gh.header.parse(encrypted_file)
    # an X25519 header packet starts with the 4 bytes of the method, then the sender's public key
    assert packet[4:36] == crypt4gh.keys.get_public_key(sender_public_key_path)
