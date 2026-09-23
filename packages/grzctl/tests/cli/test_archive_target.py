"""Tests for ``ArchiveTarget``'s public key: inline text or a path, exactly one, and
``public_key_file`` giving a file path either way.
"""

from pathlib import Path

import pytest
from grz_common.models.s3 import S3Options
from grzctl.models.config import ArchiveTarget
from pydantic import ValidationError


def _archive_target(**public_key_kwargs) -> ArchiveTarget:
    return ArchiveTarget(s3=S3Options(bucket="consented"), **public_key_kwargs)


def test_neither_public_key_nor_public_key_path_fails():
    with pytest.raises(ValidationError, match="Either public_key or public_key_path must be set"):
        _archive_target()


def test_both_public_key_and_public_key_path_fails(crypt4gh_public_key: str):
    with pytest.raises(ValidationError, match="Only one of public_key or public_key_path must be set"):
        _archive_target(public_key=crypt4gh_public_key, public_key_path="/dev/null")


def test_malformed_public_key_fails():
    with pytest.raises(ValidationError, match="Invalid public key format"):
        _archive_target(public_key="not a crypt4gh key")


def test_public_key_file_gives_the_path_unchanged():
    target = _archive_target(public_key_path="/dev/null")

    with target.public_key_file() as path:
        assert path == "/dev/null"


def test_public_key_file_writes_the_inline_key_to_a_temporary_file(crypt4gh_public_key: str):
    target = _archive_target(public_key=crypt4gh_public_key)

    with target.public_key_file() as path:
        written_path = Path(path)
        assert written_path.read_text() == crypt4gh_public_key

    assert not written_path.exists(), "the temporary file must be cleaned up once the caller is done with it"
