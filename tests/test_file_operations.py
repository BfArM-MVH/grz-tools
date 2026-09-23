"""Tests for the file_operations module."""

import errno
import os
from pathlib import Path

import grz_common.exceptions as grzexc
import pytest
from grz_common.utils.checksums import calculate_sha256
from grz_common.utils.crypt import Crypt4GH
from grz_common.utils.io import TqdmIOWrapper
from grz_common.utils.paths import is_relative_subdirectory


def test_calculate_sha256(temp_small_file_path: str, temp_small_file_sha256sum):
    sha256 = calculate_sha256(temp_small_file_path)
    assert isinstance(sha256, str)
    assert len(sha256) == 64  # sha256 hash is 64 characters long
    assert sha256 == temp_small_file_sha256sum


def test_prepare_c4gh_keys(crypt4gh_grz_public_key_file_path: str):
    keys = Crypt4GH.prepare_c4gh_keys(crypt4gh_grz_public_key_file_path)
    # single key in tuple
    assert len(keys) == 1
    # key method is set to 0
    assert keys[0][0] == 0
    # private key is generated
    assert len(keys[0][1]) == 32


@pytest.mark.parametrize(
    "relative_path, root_directory, expected",
    [
        # Valid subdirectory paths
        ("root/directory/subdir/file.bed", "root/directory", True),
        ("root/directory/subdir", "root/directory", True),
        ("root/directory/another_subdir/file.bed", "root/directory", True),
        # Target path is exactly the root directory
        ("root/directory", "root/directory", True),
        # Trying to escape root
        ("root/directory/../file_outside.bed", "root/directory", False),
        ("root/directory/../../outside/file.bed", "root/directory", False),
        # Same as root with different formatting
        ("root/directory/.", "root/directory", True),
        ("root/directory/./subdir", "root/directory", True),
        # Completely different path
        ("/some/other/directory/file.bed", "/home/user/projects/root", False),
        ("other/directory", "root/directory", False),
    ],
)
def test_is_relative_subdirectory(relative_path, root_directory, expected):
    """
    Test the is_relative_subdirectory() function with various cases.
    """
    result = is_relative_subdirectory(Path(relative_path), Path(root_directory))
    assert result == expected


def test_crypt4gh_encrypt_file(
    temp_small_file_path: str,
    crypt4gh_grz_public_keys,
    crypt4gh_grz_private_key_file_path,
    tmp_path_factory,
):
    tmp_dir = tmp_path_factory.mktemp("crypt4gh")

    tmp_encrypted_file = tmp_dir / "temp_file.c4gh"
    tmp_decrypted_file = tmp_dir / "temp_file"

    Crypt4GH.encrypt_file(temp_small_file_path, tmp_encrypted_file, crypt4gh_grz_public_keys)

    private_key = Crypt4GH.retrieve_private_key(crypt4gh_grz_private_key_file_path)

    Crypt4GH.decrypt_file(tmp_encrypted_file, tmp_decrypted_file, private_key=private_key)

    import filecmp

    assert filecmp.cmp(temp_small_file_path, tmp_decrypted_file)


def test_crypt4gh_decrypt_file_reports_a_changed_byte_as_a_decryption_error(
    temp_small_file_path: str,
    crypt4gh_grz_public_keys,
    crypt4gh_grz_private_key_file_path,
    tmp_path,
):
    """The file is at fault, so the error blames the submission and not the setup."""
    tmp_encrypted_file = tmp_path / "temp_file.c4gh"
    Crypt4GH.encrypt_file(temp_small_file_path, tmp_encrypted_file, crypt4gh_grz_public_keys)
    encrypted = bytearray(tmp_encrypted_file.read_bytes())
    # the last byte belongs to the MAC of the last segment
    encrypted[-1] ^= 0xFF
    tmp_encrypted_file.write_bytes(bytes(encrypted))
    private_key = Crypt4GH.retrieve_private_key(crypt4gh_grz_private_key_file_path)

    with pytest.raises(grzexc.DecryptionError):
        Crypt4GH.decrypt_file(tmp_encrypted_file, tmp_path / "temp_file", private_key=private_key)


@pytest.fixture
def encrypted_random_file(crypt4gh_grz_public_keys, tmp_path) -> Path:
    """A file of several crypt4gh segments."""
    plain_file = tmp_path / "random"
    plain_file.write_bytes(os.urandom(300 * 1024))
    encrypted_file = tmp_path / "random.c4gh"
    Crypt4GH.encrypt_file(plain_file, encrypted_file, crypt4gh_grz_public_keys)
    return encrypted_file


@pytest.mark.skipif(not Path("/dev/full").exists(), reason="needs /dev/full")
def test_crypt4gh_decrypt_file_raises_the_error_of_a_full_disk(
    encrypted_random_file, crypt4gh_grz_private_key_file_path
):
    """A full disk does not look like a decrypted file, although crypt4gh swallows the OSError of the failed write."""
    private_key = Crypt4GH.retrieve_private_key(crypt4gh_grz_private_key_file_path)

    with pytest.raises(OSError) as excinfo:
        Crypt4GH.decrypt_file(encrypted_random_file, Path("/dev/full"), private_key=private_key)

    assert excinfo.value.errno == errno.ENOSPC


def test_crypt4gh_decrypt_file_raises_the_error_of_a_failed_read(
    encrypted_random_file, crypt4gh_grz_private_key_file_path, tmp_path, monkeypatch
):
    """A failed read does not look like a decrypted file either."""
    private_key = Crypt4GH.retrieve_private_key(crypt4gh_grz_private_key_file_path)

    def fail(*_args):
        raise OSError(errno.EIO, "Input/output error")

    monkeypatch.setattr(TqdmIOWrapper, "readinto", fail)

    with pytest.raises(OSError) as excinfo:
        Crypt4GH.decrypt_file(encrypted_random_file, tmp_path / "random", private_key=private_key)

    assert excinfo.value.errno == errno.EIO
