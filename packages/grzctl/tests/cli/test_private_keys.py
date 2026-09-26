"""Tests for the crypt4gh private keys in the grzctl config: each is given inline as ``<name>`` or as a
path in ``<name>_path``, with an optional ``<name>_passphrase``, and loads in memory.
"""

from pathlib import Path
from unittest.mock import patch

import crypt4gh.keys
import crypt4gh.keys.c4gh
import grz_common.exceptions as grzexc
import pytest
import yaml
from grz_common.models.s3 import S3Options
from grzctl.models.config import ArchivesConfig, ArchiveTarget, GrzctlConfig, InboxConfig
from pydantic import ValidationError

from .conftest import PRUEFBERICHT

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


def _archive_target(bucket: str, unread_file: str) -> ArchiveTarget:
    return ArchiveTarget(s3=S3Options(bucket=bucket), public_key_path=unread_file)


def test_both_archive_private_key_and_private_key_path_fails(key_path: Path, unread_file: str):
    with pytest.raises(ValidationError, match="Only one of private_key or private_key_path must be set"):
        ArchiveTarget(
            s3=S3Options(bucket="consented"),
            public_key_path=unread_file,
            private_key=key_path.read_text(),
            private_key_path=str(key_path),
        )


def test_inline_archive_private_key_is_named_by_its_config_location_in_errors(no_prompt, unread_file: str):
    archives = ArchivesConfig(
        consented=_archive_target("consented", unread_file),
        non_consented=ArchiveTarget(
            s3=S3Options(bucket="non_consented"), public_key_path=unread_file, private_key="not a key"
        ),
    )

    with pytest.raises(
        grzexc.ConfigurationError, match=r"Secret key archives\.non_consented\.private_key cannot be read"
    ) as exc_info:
        archives.load_private_key("non_consented")

    assert "not a key" not in str(exc_info.value)


def _inbox(**private_key_fields) -> InboxConfig:
    return InboxConfig(**private_key_fields)


def test_neither_inbox_private_key_nor_private_key_path_fails():
    with pytest.raises(ValidationError, match="Either private_key or private_key_path must be set"):
        _inbox()


def test_both_inbox_private_key_and_private_key_path_fails(key_path: Path):
    with pytest.raises(ValidationError, match="Only one of private_key or private_key_path must be set"):
        _inbox(private_key=key_path.read_text(), private_key_path=str(key_path))


def _other_sections(tmp_path: Path, unread_file: str) -> dict:
    """The sections of a config other than ``leistungserbringer``."""
    return {
        "archives": {
            name: {"s3": {"bucket": name}, "public_key_path": unread_file} for name in ("consented", "non_consented")
        },
        "db": {"database_url": f"sqlite:///{tmp_path / 'unused.sqlite'}", "author": {"name": "test"}},
        "pruefbericht": PRUEFBERICHT,
        "identifiers": {"grz": "GRZK00007"},
    }


def _grzctl_config(tmp_path: Path, unread_file: str, leistungserbringer: dict) -> GrzctlConfig:
    """A config with the given inboxes."""
    return GrzctlConfig.from_configuration(
        {"leistungserbringer": leistungserbringer, **_other_sections(tmp_path, unread_file)}
    )


def _config_with_one_inbox(tmp_path: Path, unread_file: str, **private_key_fields) -> GrzctlConfig:
    """A config whose submitter 260914050 has one inbox, with the given private key fields."""
    return _grzctl_config(tmp_path, unread_file, {"260914050": {"inbox_buckets": {"inbox": private_key_fields}}})


def test_a_config_error_does_not_show_the_passphrase(tmp_path: Path, key_path: Path, unread_file: str):
    """Pydantic shows the start and the end of the raw input in errors, and SecretStr does not mask it there."""
    inbox = {
        "private_key": key_path.read_text(),
        "private_key_path": str(key_path),
        "private_key_passphrase": PASSPHRASE,
    }

    with pytest.raises(ValidationError, match="Only one of private_key or private_key_path must be set") as exc_info:
        _grzctl_config(tmp_path, unread_file, {"260914050": {"inbox_buckets": {"inbox": inbox}}})

    assert PASSPHRASE not in str(exc_info.value)


def test_inline_private_key_loads_in_memory(
    tmp_path: Path, key_path: Path, expected_key: bytes, no_prompt, unread_file: str
):
    config = _config_with_one_inbox(
        tmp_path, unread_file, private_key=key_path.read_text(), private_key_passphrase=PASSPHRASE
    )

    with (
        patch("builtins.open", side_effect=AssertionError("no file may be opened")),
        patch("tempfile.NamedTemporaryFile", side_effect=AssertionError("no temporary file may be written")),
        patch("tempfile.mkstemp", side_effect=AssertionError("no temporary file may be written")),
    ):
        loaded = config.inbox_target("260914050", "inbox").load_private_key()

    assert loaded.private_bytes_raw() == expected_key


def test_private_key_path_loads_with_its_passphrase(
    tmp_path: Path, key_path: Path, expected_key: bytes, no_prompt, unread_file: str
):
    config = _config_with_one_inbox(
        tmp_path, unread_file, private_key_path=str(key_path), private_key_passphrase=PASSPHRASE
    )

    assert config.inbox_target("260914050", "inbox").load_private_key().private_bytes_raw() == expected_key


@pytest.fixture
def key_paths(tmp_path: Path) -> dict[str, Path]:
    """Two crypt4gh private keys without passphrase, by name."""
    paths = {}
    for name in ("first", "second"):
        paths[name] = tmp_path / f"{name}.sec"
        crypt4gh.keys.c4gh.generate(paths[name], tmp_path / f"{name}.pub", None, comment=None)
    return paths


def test_inbox_target_loads_the_key_of_its_own_inbox(
    tmp_path: Path, key_paths: dict[str, Path], no_prompt, unread_file: str
):
    """Each inbox names its own key, so two inboxes of one submitter may use different keys."""
    config = _grzctl_config(
        tmp_path,
        unread_file,
        {
            "260914050": {
                "inbox_buckets": {
                    "inbox-a": {"private_key": key_paths["first"].read_text()},
                    "inbox-b": {"private_key_path": str(key_paths["second"])},
                }
            }
        },
    )

    loaded = config.inbox_target("260914050", "inbox-b").load_private_key()

    assert loaded.private_bytes_raw() == crypt4gh.keys.get_private_key(key_paths["second"], None)


def test_yaml_anchors_share_one_key_between_two_inboxes(
    tmp_path: Path, key_path: Path, expected_key: bytes, no_prompt, unread_file: str
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
        f"    public_key_path: {unread_file}\n"
        "  non_consented:\n"
        "    s3: {bucket: non_consented}\n"
        f"    public_key_path: {unread_file}\n"
        "db:\n"
        f"  database_url: sqlite:///{tmp_path / 'unused.sqlite'}\n"
        "  author: {name: test}\n"
        "pruefbericht:\n"
        "  authorization_url: https://auth.example.org\n"
        "  client_id: example-client\n"
        "  client_secret: example-secret\n"
        "  api_base_url: https://api.example.org\n"
        "identifiers: {grz: GRZK00007}\n"
    )

    config = GrzctlConfig.from_path(config_path)

    for inbox_name in ("inbox-a", "inbox-b"):
        loaded = config.inbox_target("260914050", inbox_name).load_private_key()
        assert loaded.private_bytes_raw() == expected_key, inbox_name


INBOX_DEFAULTS_LAYOUT = """\
inbox_defaults: &inbox_defaults # grzctl ignores this section; the inboxes merge it with <<
  endpoint_url: https://s3.example.org
  private_key_path: {private_key_path}
leistungserbringer:
  "123456789": # LE ID
    alias: "FOO" # optional
    inbox_buckets:
      main: # the inbox name, which --inbox takes
        <<: *inbox_defaults
        bucket: le-123456789 # optional, defaults to the inbox name
  "000000000":
    inbox_buckets:
      main:
        <<: *inbox_defaults
        bucket: le-000000000
"""
"""The ``inbox_defaults`` layout of ``crypt4gh-keys.md``, with a placeholder for the key file."""


def test_inboxes_merge_the_inbox_defaults(tmp_path: Path, key_path: Path, unread_file: str):
    """Each inbox merges the top-level ``inbox_defaults`` with ``<<``, as ``crypt4gh-keys.md`` shows.

    ``grzctl --config`` loads the file with :meth:`GrzctlConfig.from_path`.
    """
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        INBOX_DEFAULTS_LAYOUT.format(private_key_path=key_path) + yaml.safe_dump(_other_sections(tmp_path, unread_file))
    )

    config = GrzctlConfig.from_path(config_path)

    for submitter_id in ("123456789", "000000000"):
        target = config.inbox_target(submitter_id, "main")
        assert str(target.s3.endpoint_url) == "https://s3.example.org/", submitter_id
        assert target.private_key_path == key_path, submitter_id
        assert target.s3.bucket == f"le-{submitter_id}", submitter_id


def test_a_config_with_the_removed_signing_key_fields_still_loads(tmp_path: Path, unread_file: str):
    """The config ignores ``archives.signing_key``, ``signing_key_path`` and ``signing_key_passphrase``.

    ``grzctl encrypt`` signs with the key of the submission's inbox instead.
    Not even a missing file in ``signing_key_path`` fails the config.
    """
    data = _other_sections(tmp_path, unread_file)
    data["leistungserbringer"] = {"260914050": {"inbox_buckets": {"inbox": {"private_key_path": unread_file}}}}
    data["archives"]["signing_key"] = "not a key"
    data["archives"]["signing_key_path"] = str(tmp_path / "missing.sec")
    data["archives"]["signing_key_passphrase"] = PASSPHRASE
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(data))

    config = GrzctlConfig.from_path(config_path)

    assert "signing_key_path" not in config.archives.model_dump()


@pytest.mark.parametrize(
    ("section", "field"),
    [
        (("leistungserbringer", "260914050", "inbox_buckets", "inbox"), "private_key_path"),
        (("archives", "consented"), "public_key_path"),
    ],
)
def test_a_missing_key_file_fails_loading_the_config(
    tmp_path: Path, unread_file: str, section: tuple[str, ...], field: str
):
    """A key path must name an existing file, so a missing file fails the config for every command."""
    inbox = {"private_key_path": unread_file}
    config = _grzctl_config(tmp_path, unread_file, {"260914050": {"inbox_buckets": {"inbox": inbox}}})
    data = config.model_dump(mode="json", exclude_none=True)
    fields = data
    for key in section:
        fields = fields[key]
    fields[field] = str(tmp_path / "missing.sec")

    with pytest.raises(ValidationError, match="Path does not point to a file"):
        GrzctlConfig.from_configuration(data)


def test_a_key_path_expands_the_home_directory(
    tmp_path: Path, unread_file: str, key_path: Path, expected_key: bytes, monkeypatch, no_prompt
):
    monkeypatch.setenv("HOME", str(key_path.parent))

    config = _config_with_one_inbox(
        tmp_path, unread_file, private_key_path=f"~/{key_path.name}", private_key_passphrase=PASSPHRASE
    )

    assert config.leistungserbringer["260914050"].inbox_buckets["inbox"].private_key_path == key_path
    assert config.inbox_target("260914050", "inbox").load_private_key().private_bytes_raw() == expected_key


def test_inbox_target_fails_for_a_submitter_missing_from_the_config(tmp_path: Path, unread_file: str):
    config = _grzctl_config(
        tmp_path, unread_file, {"111111111": {"inbox_buckets": {"inbox": {"private_key_path": unread_file}}}}
    )

    with pytest.raises(grzexc.ConfigurationError, match=r"Submitter '260914050' not found\. Available: '111111111'"):
        config.inbox_target("260914050", "inbox")


def test_inbox_target_fails_for_an_inbox_missing_from_the_config(tmp_path: Path, unread_file: str):
    config = _config_with_one_inbox(tmp_path, unread_file, private_key_path=unread_file)

    with pytest.raises(
        grzexc.ConfigurationError, match=r"Inbox 'other' not configured for submitter '260914050'\. Available: inbox"
    ):
        config.inbox_target("260914050", "other")


def test_inline_inbox_key_asks_for_its_passphrase_by_its_config_location(
    tmp_path: Path, key_path: Path, expected_key: bytes, monkeypatch, unread_file: str
):
    monkeypatch.delenv("C4GH_PASSPHRASE", raising=False)
    prompts = []

    def _getpass(prompt):
        prompts.append(prompt)
        return PASSPHRASE

    monkeypatch.setattr("grz_common.utils.crypt.getpass", _getpass)
    config = _config_with_one_inbox(tmp_path, unread_file, private_key=key_path.read_text())

    loaded = config.inbox_target("260914050", "inbox").load_private_key()

    assert loaded.private_bytes_raw() == expected_key
    assert prompts == ["Passphrase for leistungserbringer.260914050.inbox_buckets.inbox.private_key: "]


def test_inline_inbox_key_is_named_by_its_config_location_in_errors(tmp_path: Path, no_prompt, unread_file: str):
    config = _config_with_one_inbox(tmp_path, unread_file, private_key="not a key")
    inbox = config.inbox_target("260914050", "inbox")

    with pytest.raises(
        grzexc.ConfigurationError,
        match=r"Secret key leistungserbringer\.260914050\.inbox_buckets\.inbox\.private_key cannot be read",
    ) as exc_info:
        inbox.load_private_key()

    assert "not a key" not in str(exc_info.value)
