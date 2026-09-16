"""Tests for reading the ``known_public_keys`` file of ``grzctl db``."""

from pathlib import Path

import pytest
from cryptography.hazmat.primitives.asymmetric.ed25519 import Ed25519PrivateKey
from cryptography.hazmat.primitives.serialization import Encoding, PublicFormat
from grz_db.errors import DatabaseConfigurationError
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
