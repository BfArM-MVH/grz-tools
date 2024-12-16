"""Tests for the file_operations module."""

from pathlib import Path

import pytest

from grz_cli.file_operations import (
    Crypt4GH,
    calculate_md5,
    calculate_sha256,
    is_relative_subdirectory,
)


def test_calculate_sha256(temp_small_file_path: str, temp_small_file_sha256sum):
    sha256 = calculate_sha256(temp_small_file_path)
    assert isinstance(sha256, str)
    assert len(sha256) == 64  # sha256 hash is 64 characters long
    assert sha256 == temp_small_file_sha256sum


def test_calculate_md5(temp_small_file_path: str, temp_small_file_md5sum):
    md5 = calculate_md5(temp_small_file_path)
    assert isinstance(md5, str)
    assert len(md5) == 32  # MD5 hash is 32 characters long
    assert md5 == temp_small_file_md5sum


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
