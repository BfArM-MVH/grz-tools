"""Tests for the crypt4gh private keys in the grzctl config: each is given inline as ``<name>`` or as a
path in ``<name>_path``, with an optional ``<name>_passphrase``, and loads in memory.
"""

from pathlib import Path
from unittest.mock import patch

import crypt4gh.keys
import crypt4gh.keys.c4gh
import grz_common.exceptions as grzexc
import pytest
from grz_common.models.s3 import S3Options
from grz_common.utils.crypt import Crypt4GH
from grzctl.models.config import ArchivesConfig, ArchiveTarget, GrzctlConfig, InboxConfig
from pydantic import ValidationError

PASSPHRASE = "grz-key-passphrase"


@pytest.fixture
def key_path(tmp_path: Path) -> Path:
    """A crypt4gh private key, encrypted with ``PASSPHRASE``."""
    private_key_path = tmp_path / "grz.sec"
    crypt4gh.keys.c4gh.generate(private_key_path, tmp_path / "grz.pub", PASSPHRASE.encode(), comment=None)
    return private_key_path


@pytest.fixture
def expected_key(key_path: Path) -> bytes:
    return crypt4gh.keys.get_private_key(key_path, lambda: PASSPHRASE)


@pytest.fixture
def no_prompt(monkeypatch):
    """Fail the test if the passphrase prompt opens, and set a wrong ``C4GH_PASSPHRASE``.

    The configured passphrase comes first, so the wrong one in the environment must not matter.
    """
    monkeypatch.setenv("C4GH_PASSPHRASE", "wrong-passphrase")

    def _fail(*args, **kwargs):
        raise AssertionError("the passphrase prompt must not open")

    monkeypatch.setattr("grz_common.utils.crypt.getpass", _fail)


def _archive_target(bucket: str) -> ArchiveTarget:
    return ArchiveTarget(s3=S3Options(bucket=bucket), public_key_path="/dev/null")


def _archives(**signing_key_fields) -> ArchivesConfig:
    return ArchivesConfig(
        consented=_archive_target("consented"),
        non_consented=_archive_target("non_consented"),
        **signing_key_fields,
    )


def test_neither_signing_key_nor_signing_key_path_fails():
    with pytest.raises(ValidationError, match="Either signing_key or signing_key_path must be set"):
        _archives()


def test_both_signing_key_and_signing_key_path_fails(key_path: Path):
    with pytest.raises(ValidationError, match="Only one of signing_key or signing_key_path must be set"):
        _archives(signing_key=key_path.read_text(), signing_key_path=str(key_path))


def test_signing_key_path_loads_with_its_passphrase(key_path: Path, expected_key: bytes, no_prompt):
    archives = _archives(signing_key_path=str(key_path), signing_key_passphrase=PASSPHRASE)

    assert archives.load_signing_key() == expected_key


def test_inline_signing_key_loads_in_memory(key_path: Path, expected_key: bytes, no_prompt):
    archives = _archives(signing_key=key_path.read_text(), signing_key_passphrase=PASSPHRASE)

    with (
        patch("builtins.open", side_effect=AssertionError("no file may be opened")),
        patch("tempfile.NamedTemporaryFile", side_effect=AssertionError("no temporary file may be written")),
        patch("tempfile.mkstemp", side_effect=AssertionError("no temporary file may be written")),
    ):
        loaded = archives.load_signing_key()

    assert loaded == expected_key


def test_inline_signing_key_is_named_by_its_config_location_in_errors(no_prompt):
    archives = _archives(signing_key="not a key")

    with pytest.raises(grzexc.ConfigurationError, match=r"Secret key archives\.signing_key cannot be read") as exc_info:
        archives.load_signing_key()

    assert "not a key" not in str(exc_info.value)


def _inbox(**private_key_fields) -> InboxConfig:
    return InboxConfig(**private_key_fields)


def _archive_target_with_private_key(**private_key_fields) -> ArchiveTarget:
    return ArchiveTarget(s3=S3Options(bucket="consented"), public_key_path="/dev/null", **private_key_fields)


def test_neither_inbox_private_key_nor_private_key_path_fails():
    with pytest.raises(ValidationError, match="Either private_key or private_key_path must be set"):
        _inbox()


def test_both_inbox_private_key_and_private_key_path_fails(key_path: Path):
    with pytest.raises(ValidationError, match="Only one of private_key or private_key_path must be set"):
        _inbox(private_key=key_path.read_text(), private_key_path=str(key_path))


def test_archive_private_key_is_optional():
    target = _archive_target_with_private_key()

    assert target.private_key is None
    assert target.private_key_path is None


def test_both_archive_private_key_and_private_key_path_fails(key_path: Path):
    with pytest.raises(ValidationError, match="Only one of private_key or private_key_path must be set"):
        _archive_target_with_private_key(private_key=key_path.read_text(), private_key_path=str(key_path))


def _grzctl_config(tmp_path: Path, leistungserbringer: dict, **archive_private_keys: dict) -> GrzctlConfig:
    """A config with the given inboxes, and the archive private keys given by archive name."""
    archives = {
        name: {"s3": {"bucket": name}, "public_key_path": "/dev/null", **archive_private_keys.get(name, {})}
        for name in ("consented", "non_consented")
    }
    return GrzctlConfig.from_configuration(
        {
            "leistungserbringer": leistungserbringer,
            "archives": {**archives, "signing_key_path": "/dev/null"},
            "db": {"database_url": f"sqlite:///{tmp_path / 'unused.sqlite'}", "author": {"name": "test"}},
            "pruefbericht": {},
            "identifiers": {"grz": "GRZK00007"},
        }
    )


def _config_with_one_private_key(tmp_path: Path, holder: str, **private_key_fields) -> GrzctlConfig:
    """A config whose only key for submitter 260914050 is in its inbox, or in the consented archive."""
    if holder == "inbox":
        return _grzctl_config(tmp_path, {"260914050": {"inbox_buckets": {"inbox": private_key_fields}}})
    other_le = {"111111111": {"inbox_buckets": {"other": {"private_key_path": "/dev/null"}}}}
    return _grzctl_config(tmp_path, other_le, consented=private_key_fields)


@pytest.mark.parametrize("holder", ["inbox", "archive"])
def test_inline_private_key_loads_in_memory(
    holder: str, tmp_path: Path, key_path: Path, expected_key: bytes, no_prompt
):
    config = _config_with_one_private_key(
        tmp_path, holder, private_key=key_path.read_text(), private_key_passphrase=PASSPHRASE
    )

    with (
        patch("builtins.open", side_effect=AssertionError("no file may be opened")),
        patch("tempfile.NamedTemporaryFile", side_effect=AssertionError("no temporary file may be written")),
        patch("tempfile.mkstemp", side_effect=AssertionError("no temporary file may be written")),
    ):
        keys = list(config.iter_decryption_keys("260914050"))

    assert [key for _, key in keys] == [expected_key]


@pytest.mark.parametrize("holder", ["inbox", "archive"])
def test_private_key_path_loads_with_its_passphrase(
    holder: str, tmp_path: Path, key_path: Path, expected_key: bytes, no_prompt
):
    config = _config_with_one_private_key(
        tmp_path, holder, private_key_path=str(key_path), private_key_passphrase=PASSPHRASE
    )

    assert [key for _, key in config.iter_decryption_keys("260914050")] == [expected_key]


@pytest.fixture
def key_paths(tmp_path: Path) -> dict[str, Path]:
    """Four crypt4gh private keys without passphrase, by name."""
    paths = {}
    for name in ("first", "second", "third", "fourth"):
        paths[name] = tmp_path / f"{name}.sec"
        crypt4gh.keys.c4gh.generate(paths[name], tmp_path / f"{name}.pub", None, comment=None)
    return paths


def _load(path: Path) -> bytes:
    return crypt4gh.keys.get_private_key(path, None)


def test_iter_decryption_keys_gives_the_submitter_inboxes_first_then_the_archives(
    tmp_path: Path, key_paths: dict[str, Path], no_prompt
):
    config = _grzctl_config(
        tmp_path,
        {
            "111111111": {"inbox_buckets": {"other": {"private_key_path": "/nonexistent/other.sec"}}},
            "260914050": {
                "inbox_buckets": {
                    "inbox-a": {"private_key": key_paths["first"].read_text()},
                    "inbox-b": {"private_key_path": str(key_paths["second"])},
                }
            },
        },
        consented={"private_key_path": str(key_paths["third"])},
        non_consented={"private_key": key_paths["fourth"].read_text()},
    )

    keys = list(config.iter_decryption_keys("260914050"))

    assert keys == [
        ("leistungserbringer.260914050.inbox_buckets.inbox-a.private_key", _load(key_paths["first"])),
        ("leistungserbringer.260914050.inbox_buckets.inbox-b.private_key_path", _load(key_paths["second"])),
        ("archives.consented.private_key_path", _load(key_paths["third"])),
        ("archives.non_consented.private_key", _load(key_paths["fourth"])),
    ]


def test_iter_decryption_keys_loads_a_key_only_when_asked_for(tmp_path: Path, key_paths: dict[str, Path], no_prompt):
    config = _grzctl_config(
        tmp_path,
        {"260914050": {"inbox_buckets": {"inbox": {"private_key_path": str(key_paths["first"])}}}},
        consented={"private_key_path": "/nonexistent/archive.sec"},
    )

    keys = config.iter_decryption_keys("260914050")

    assert next(keys) == (
        "leistungserbringer.260914050.inbox_buckets.inbox.private_key_path",
        _load(key_paths["first"]),
    )
    with pytest.raises(grzexc.ConfigurationError, match=r"archives\.consented\.private_key_path: Secret key not found"):
        next(keys)


def test_iter_decryption_keys_skips_the_inboxes_of_a_submitter_missing_from_the_config(
    tmp_path: Path, key_paths: dict[str, Path], no_prompt
):
    config = _grzctl_config(
        tmp_path,
        {"111111111": {"inbox_buckets": {"other": {"private_key_path": "/nonexistent/other.sec"}}}},
        consented={"private_key_path": str(key_paths["first"])},
    )

    keys = list(config.iter_decryption_keys("260914050"))

    assert keys == [("archives.consented.private_key_path", _load(key_paths["first"]))]


def test_yaml_anchors_share_one_key_between_two_inboxes_and_the_signing_key(
    tmp_path: Path, key_path: Path, expected_key: bytes, no_prompt
):
    """A GRZ with one key pair for all inboxes writes the key once and refers to it with YAML aliases."""
    key_block = "".join(f"          {line}\n" for line in key_path.read_text().splitlines())
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        "leistungserbringer:\n"
        "  '260914050':\n"
        "    inbox_buckets:\n"
        "      inbox-a:\n"
        "        private_key: &grz_key |\n"
        f"{key_block}"
        f"        private_key_passphrase: &grz_key_passphrase {PASSPHRASE}\n"
        "      inbox-b:\n"
        "        private_key: *grz_key\n"
        "        private_key_passphrase: *grz_key_passphrase\n"
        "archives:\n"
        "  consented:\n"
        "    s3: {bucket: consented}\n"
        "    public_key_path: /dev/null\n"
        "  non_consented:\n"
        "    s3: {bucket: non_consented}\n"
        "    public_key_path: /dev/null\n"
        "  signing_key: *grz_key\n"
        "  signing_key_passphrase: *grz_key_passphrase\n"
        "db:\n"
        f"  database_url: sqlite:///{tmp_path / 'unused.sqlite'}\n"
        "  author: {name: test}\n"
        "pruefbericht: {}\n"
        "identifiers: {grz: GRZK00007}\n"
    )

    config = GrzctlConfig.from_path(config_path)

    with patch.object(Crypt4GH, "load_private_key", wraps=Crypt4GH.load_private_key) as load_private_key:
        keys = list(config.iter_decryption_keys("260914050"))

    assert keys == [
        (
            "leistungserbringer.260914050.inbox_buckets.inbox-a.private_key, "
            "leistungserbringer.260914050.inbox_buckets.inbox-b.private_key",
            expected_key,
        )
    ]
    assert load_private_key.call_count == 1, "the inline key that both inboxes share loads once"
    assert config.archives.load_signing_key() == expected_key


def test_two_inboxes_sharing_a_key_path_through_a_yaml_anchor_prompt_once(
    tmp_path: Path, key_path: Path, expected_key: bytes, monkeypatch
):
    """Without a configured passphrase, the shared key asks for its passphrase once, not once per inbox."""
    monkeypatch.delenv("C4GH_PASSPHRASE", raising=False)
    prompts = []

    def _getpass(prompt):
        prompts.append(prompt)
        return PASSPHRASE

    monkeypatch.setattr("grz_common.utils.crypt.getpass", _getpass)
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        "leistungserbringer:\n"
        "  '260914050':\n"
        "    inbox_buckets:\n"
        "      inbox-a:\n"
        f"        private_key_path: &grz_key {key_path}\n"
        "      inbox-b:\n"
        "        private_key_path: *grz_key\n"
        "archives:\n"
        "  consented:\n"
        "    s3: {bucket: consented}\n"
        "    public_key_path: /dev/null\n"
        "  non_consented:\n"
        "    s3: {bucket: non_consented}\n"
        "    public_key_path: /dev/null\n"
        "  signing_key_path: *grz_key\n"
        "db:\n"
        f"  database_url: sqlite:///{tmp_path / 'unused.sqlite'}\n"
        "  author: {name: test}\n"
        "pruefbericht: {}\n"
        "identifiers: {grz: GRZK00007}\n"
    )

    keys = list(GrzctlConfig.from_path(config_path).iter_decryption_keys("260914050"))

    assert keys == [
        (
            "leistungserbringer.260914050.inbox_buckets.inbox-a.private_key_path, "
            "leistungserbringer.260914050.inbox_buckets.inbox-b.private_key_path",
            expected_key,
        )
    ]
    assert prompts == [f"Passphrase for {key_path}: "]
