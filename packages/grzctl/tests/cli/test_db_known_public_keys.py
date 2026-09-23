"""Tests for the known public keys of ``grzctl db`` (inline or from a file), and for the signature
checks that use them.
"""

from pathlib import Path

import click.testing
import grzctl.cli
import pytest
from cryptography.hazmat.primitives.asymmetric.ed25519 import Ed25519PrivateKey
from cryptography.hazmat.primitives.serialization import Encoding, PublicFormat
from grz_db.errors import DatabaseConfigurationError
from grzctl.commands.db import SignatureStatus, _verify_signature
from grzctl.commands.db.cli import _parse_known_public_keys, _read_known_public_keys_file
from grzctl.models.config import GrzctlConfig
from grzctl.models.db import Author, DbModel
from pydantic import ValidationError


def _openssh_public_key() -> str:
    """A fresh ed25519 public key as ``<format> <key>``, without a comment."""
    public_key = Ed25519PrivateKey.generate().public_key()
    return public_key.public_bytes(encoding=Encoding.OpenSSH, format=PublicFormat.OpenSSH).decode()


def _write(tmp_path: Path, content: str) -> Path:
    path = tmp_path / "known_public_keys"
    path.write_text(content)
    return path


def test_file_skips_blank_lines_and_comment_lines(tmp_path: Path) -> None:
    path = _write(tmp_path, f"# data stewards\n\n{_openssh_public_key()} alice\n\n{_openssh_public_key()} bob\n")

    assert _read_known_public_keys_file(path).keys() == {"alice", "bob"}


def test_file_names_the_line_of_an_error(tmp_path: Path) -> None:
    path = _write(tmp_path, f"# data stewards\n{_openssh_public_key()}\n")

    with pytest.raises(DatabaseConfigurationError, match=r"known_public_keys:2: expected"):
        _read_known_public_keys_file(path)


def test_keeps_every_key_that_shares_a_comment_in_order() -> None:
    """A rotated key keeps its owner's name, and the owner's older signatures must still verify."""
    first = _openssh_public_key()
    second = _openssh_public_key()
    entries = [("db.known_public_keys[0]", f"{first} alice"), ("db.known_public_keys[1]", f"{second} alice")]

    keys = _parse_known_public_keys(entries)["alice"]

    assert [key.public_bytes(Encoding.OpenSSH, PublicFormat.OpenSSH).decode() for key in keys] == [first, second]


def test_keeps_a_comment_with_spaces_whole() -> None:
    entries = [("db.known_public_keys[0]", f"{_openssh_public_key()} Alice Example")]

    assert _parse_known_public_keys(entries).keys() == {"Alice Example"}


@pytest.mark.parametrize("entry", ["", "# data stewards"], ids=["blank", "comment"])
def test_rejects_an_entry_that_is_no_key_and_names_it(entry: str) -> None:
    """The list takes its comments from YAML, so an entry that is no key is a mistake."""
    entries = [("db.known_public_keys[0]", f"{_openssh_public_key()} alice"), ("db.known_public_keys[1]", entry)]

    with pytest.raises(DatabaseConfigurationError, match=r"db\.known_public_keys\[1\]: "):
        _parse_known_public_keys(entries)


def test_rejects_a_key_without_a_comment_and_names_it() -> None:
    with pytest.raises(DatabaseConfigurationError, match=r"db\.known_public_keys\[0\]: expected"):
        _parse_known_public_keys([("db.known_public_keys[0]", _openssh_public_key())])


def test_rejects_a_key_that_does_not_load_and_names_it() -> None:
    with pytest.raises(DatabaseConfigurationError, match=r"db\.known_public_keys\[0\]: cannot load"):
        _parse_known_public_keys([("db.known_public_keys[0]", "ssh-ed25519 not-a-key alice")])


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


def _db_model(**known_public_keys_kwargs) -> DbModel:
    return DbModel(
        database_url="sqlite:///unused.sqlite",
        author=Author(name="alice", private_key="dummy"),
        **known_public_keys_kwargs,
    )


def test_known_public_keys_rejects_a_path_with_a_migration_hint() -> None:
    """Old configs put a path into ``known_public_keys``; the error points at the new field."""
    with pytest.raises(ValidationError, match="known_public_keys_file"):
        _db_model(known_public_keys="/some/known_public_keys")


def test_known_public_keys_rejects_both_set(tmp_path: Path) -> None:
    path = _write(tmp_path, f"{_openssh_public_key()} alice\n")

    with pytest.raises(ValidationError, match="Only one of known_public_keys or known_public_keys_file"):
        _db_model(known_public_keys=[f"{_openssh_public_key()} alice"], known_public_keys_file=str(path))


def test_known_public_keys_file_expands_home(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("HOME", str(tmp_path))
    path = tmp_path / "known_public_keys"
    path.write_text(f"{_openssh_public_key()} alice\n")

    config = _db_model(known_public_keys_file="~/known_public_keys")

    assert config.known_public_keys_file == path


def _write_config(tmp_path: Path, config: GrzctlConfig) -> Path:
    config_path = tmp_path / "config.yaml"
    with open(config_path, "w") as config_file:
        config.to_yaml(config_file)
    return config_path


def test_db_group_requires_one_of_known_public_keys_or_known_public_keys_file(
    tmp_path: Path, offline_config: GrzctlConfig
) -> None:
    offline_config.db.known_public_keys_file = None
    config_path = _write_config(tmp_path, offline_config)
    cli = grzctl.cli.build_cli()

    result = click.testing.CliRunner().invoke(cli, ["--config", str(config_path), "db", "init"])

    assert isinstance(result.exception, DatabaseConfigurationError)
    assert "known_public_keys" in str(result.exception)


def test_db_group_accepts_an_inline_known_public_keys_list(tmp_path: Path, offline_config: GrzctlConfig) -> None:
    """A config with the keys inlined works just as well as one naming a file."""
    keys_file = Path(offline_config.db.known_public_keys_file)
    offline_config.db.known_public_keys_file = None
    offline_config.db.known_public_keys = [line for line in keys_file.read_text().splitlines() if line.strip()]
    config_path = _write_config(tmp_path, offline_config)
    cli = grzctl.cli.build_cli()

    result = click.testing.CliRunner().invoke(cli, ["--config", str(config_path), "db", "init"])

    assert result.exit_code == 0, result.stderr
