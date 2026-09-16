"""Tests for reading the ``known_public_keys`` file of ``grzctl db``, and for the signature checks that use it."""

from pathlib import Path

import pytest
from cryptography.hazmat.primitives.asymmetric.ed25519 import Ed25519PrivateKey
from cryptography.hazmat.primitives.serialization import Encoding, PublicFormat
from grz_db.errors import DatabaseConfigurationError
from grzctl.commands.db import SignatureStatus, _verify_signature
from grzctl.commands.db.cli import _read_known_public_keys


def _openssh_public_key() -> str:
    """A fresh ed25519 public key as ``<format> <key>``, without a comment."""
    public_key = Ed25519PrivateKey.generate().public_key()
    return public_key.public_bytes(encoding=Encoding.OpenSSH, format=PublicFormat.OpenSSH).decode()


def _write(tmp_path: Path, content: str) -> Path:
    path = tmp_path / "known_public_keys"
    path.write_text(content)
    return path


def test_skips_blank_lines_and_comment_lines(tmp_path: Path) -> None:
    path = _write(tmp_path, f"# data stewards\n\n{_openssh_public_key()} alice\n\n{_openssh_public_key()} bob\n")

    assert _read_known_public_keys(path).keys() == {"alice", "bob"}


def test_keeps_every_key_that_shares_a_comment_in_file_order(tmp_path: Path) -> None:
    """A rotated key keeps its owner's name, and the owner's older signatures must still verify."""
    first = _openssh_public_key()
    second = _openssh_public_key()
    path = _write(tmp_path, f"{first} alice\n{second} alice\n")

    keys = _read_known_public_keys(path)["alice"]

    assert [key.public_bytes(Encoding.OpenSSH, PublicFormat.OpenSSH).decode() for key in keys] == [first, second]


def test_keeps_a_comment_with_spaces_whole(tmp_path: Path) -> None:
    path = _write(tmp_path, f"{_openssh_public_key()} Alice Example\n")

    assert _read_known_public_keys(path).keys() == {"Alice Example"}


def test_rejects_a_key_without_a_comment_and_names_the_line(tmp_path: Path) -> None:
    path = _write(tmp_path, f"# data stewards\n{_openssh_public_key()}\n")

    with pytest.raises(DatabaseConfigurationError, match=r"known_public_keys:2: expected"):
        _read_known_public_keys(path)


def test_rejects_a_key_that_does_not_load_and_names_the_line(tmp_path: Path) -> None:
    path = _write(tmp_path, "ssh-ed25519 not-a-key alice\n")

    with pytest.raises(DatabaseConfigurationError, match=r"known_public_keys:1: cannot load"):
        _read_known_public_keys(path)


class _SignedBy:
    """Stands in for a signed log entry that only the given public key verifies."""

    def __init__(self, public_key):
        self.public_key = public_key

    def verify(self, public_key) -> bool:
        return public_key is self.public_key


def _public_key():
    return Ed25519PrivateKey.generate().public_key()


def test_a_name_with_several_keys_verifies_with_any_of_them() -> None:
    old = _public_key()
    new = _public_key()

    assert _verify_signature({"alice": [new, old]}, "alice", _SignedBy(old)) == (SignatureStatus.VERIFIED, None)


def test_a_name_whose_keys_all_fail_reports_failed() -> None:
    public_keys = {"alice": [_public_key(), _public_key()]}

    assert _verify_signature(public_keys, "alice", _SignedBy(_public_key())) == (SignatureStatus.FAILED, None)


def test_an_unknown_name_tries_every_key_of_every_name() -> None:
    signer = _public_key()
    public_keys = {"alice": [_public_key()], "bob": [_public_key(), signer]}

    assert _verify_signature(public_keys, "carol", _SignedBy(signer)) == (SignatureStatus.VERIFIED, "bob")
