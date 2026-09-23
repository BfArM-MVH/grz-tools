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
from grzctl.models.config import ArchivesConfig, ArchiveTarget
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
