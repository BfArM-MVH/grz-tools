"""Tests for ``ArchiveTarget``'s public key: inline text or a path, exactly one, and
``load_public_key`` loading it either way.
"""

from base64 import b64decode
from pathlib import Path
from unittest.mock import patch

import grz_common.exceptions as grzexc
import pytest
from grz_common.models.s3 import S3Options
from grzctl.models.config import ArchiveTarget
from pydantic import ValidationError


def _archive_target(**public_key_kwargs) -> ArchiveTarget:
    return ArchiveTarget(s3=S3Options(bucket="consented"), **public_key_kwargs)


def test_neither_public_key_nor_public_key_path_fails():
    with pytest.raises(ValidationError, match="Either public_key or public_key_path must be set"):
        _archive_target()


def test_both_public_key_and_public_key_path_fails(crypt4gh_public_key: str, unread_file: str):
    with pytest.raises(ValidationError, match="Only one of public_key or public_key_path must be set"):
        _archive_target(public_key=crypt4gh_public_key, public_key_path=unread_file)


def test_malformed_public_key_fails():
    with pytest.raises(ValidationError, match="Invalid public key format"):
        _archive_target(public_key="not a crypt4gh key")


def test_load_public_key_loads_an_inline_key_in_memory(crypt4gh_public_key: str):
    target = _archive_target(public_key=crypt4gh_public_key)

    with (
        patch("builtins.open", side_effect=AssertionError("no file may be opened")),
        patch("tempfile.NamedTemporaryFile", side_effect=AssertionError("no temporary file may be written")),
        patch("tempfile.mkstemp", side_effect=AssertionError("no temporary file may be written")),
    ):
        loaded = target.load_public_key()

    assert loaded.public_bytes_raw() == b64decode(crypt4gh_public_key.splitlines()[1])


def test_load_public_key_loads_a_key_file(tmp_path: Path, crypt4gh_public_key: str):
    public_key_path = tmp_path / "archive.pub"
    public_key_path.write_text(crypt4gh_public_key)
    target = _archive_target(public_key_path=str(public_key_path))

    loaded = target.load_public_key()

    assert loaded.public_bytes_raw() == b64decode(crypt4gh_public_key.splitlines()[1])


def test_load_public_key_names_the_archive_of_an_inline_key_that_does_not_load():
    """``Crypt4GHPublicKey`` only checks the markers, so a payload that is no key fails only on loading."""
    target = _archive_target(
        public_key="-----BEGIN CRYPT4GH PUBLIC KEY-----\nbm90IGEga2V5\n-----END CRYPT4GH PUBLIC KEY-----\n"
    )

    with pytest.raises(
        grzexc.ConfigurationError, match="Public key of the archive with bucket consented cannot be read"
    ):
        target.load_public_key()
