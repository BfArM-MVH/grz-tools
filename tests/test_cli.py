"""Tests for the command line interface."""

import filecmp
import os
import shutil
from pathlib import Path
from unittest import mock

import pytest
from click.testing import CliRunner
from moto import mock_aws

import grz_cli.cli
from grz_cli.file_operations import calculate_sha256


@pytest.fixture
def working_dir(tmpdir_factory: pytest.TempdirFactory):
    """Create temporary folder for the session"""
    datadir = tmpdir_factory.mktemp("submission")
    return datadir


@pytest.fixture
def working_dir_path(working_dir) -> Path:
    return Path(working_dir.strpath)


def test_validate_submission(
    working_dir_path,
    temp_config_file_path,
):
    submission_dir = Path("tests/mock_files/submissions/valid_submission")

    shutil.copytree(submission_dir / "files", working_dir_path / "files", dirs_exist_ok=True)
    shutil.copytree(submission_dir / "metadata", working_dir_path / "metadata", dirs_exist_ok=True)

    testargs = [
        "validate",
        "--submission-dir",
        str(working_dir_path),
    ]

    runner = CliRunner()
    result = runner.invoke(grz_cli.cli.cli, testargs, catch_exceptions=False)

    # check if re-validation is skipped
    result = runner.invoke(grz_cli.cli.cli, testargs, catch_exceptions=False)

    # test if command has correctly checked for:
    # - mismatched md5sums
    # - all files existing

    assert result.exit_code == 0, result.output


def test_encrypt_decrypt_submission(
    working_dir_path,
    temp_config_file_path,
    # crypt4gh_grz_private_key_file_path,
    tmpdir_factory: pytest.TempdirFactory,
):
    submission_dir = Path("tests/mock_files/submissions/valid_submission")

    shutil.copytree(submission_dir / "files", working_dir_path / "files", dirs_exist_ok=True)
    shutil.copytree(submission_dir / "metadata", working_dir_path / "metadata", dirs_exist_ok=True)

    # first, encrypt the data
    testargs = [
        "encrypt",
        "--submission-dir",
        str(working_dir_path),
        "--config-file",
        temp_config_file_path,
    ]

    runner = CliRunner()
    result = runner.invoke(grz_cli.cli.cli, testargs, catch_exceptions=False)

    assert result.exit_code == 0, result.output

    # then, decrypt the data again
    testargs = [
        "decrypt",
        "--submission-dir",
        str(working_dir_path),
        "--config-file",
        temp_config_file_path,
    ]

    runner = CliRunner()
    result = runner.invoke(grz_cli.cli.cli, testargs, catch_exceptions=False)

    assert result.exit_code == 0, result.output

    # compare if the files are equal
    for file in [
        "aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000_blood_normal.read1.fastq.gz",
        "aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000_blood_normal.read2.fastq.gz",
        "aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000_blood_tumor.read1.fastq.gz",
        "aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000_blood_tumor.read2.fastq.gz",
        "bbbbbbbb11111111bbbbbbbb11111111bbbbbbbb11111111bbbbbbbb11111111_blood_normal.read1.fastq.gz",
        "bbbbbbbb11111111bbbbbbbb11111111bbbbbbbb11111111bbbbbbbb11111111_blood_normal.read2.fastq.gz",
    ]:
        expected_checksum = calculate_sha256(submission_dir / "files" / file)
        observed_checksum = calculate_sha256(working_dir_path / "files" / file)

        assert expected_checksum == observed_checksum


def test_decrypt_submission(working_dir_path, temp_config_file_path):
    submission_dir = Path("tests/mock_files/submissions/valid_submission")

    shutil.copytree(
        submission_dir / "encrypted_files",
        working_dir_path / "encrypted_files",
        dirs_exist_ok=True,
    )
    shutil.copytree(submission_dir / "metadata", working_dir_path / "metadata", dirs_exist_ok=True)

    testargs = [
        "decrypt",
        "--submission-dir",
        str(working_dir_path),
        "--config-file",
        temp_config_file_path,
    ]
    runner = CliRunner()
    result = runner.invoke(grz_cli.cli.cli, testargs, catch_exceptions=False)

    assert result.exit_code == 0, result.output

    # compare if the files are equal
    for file in [
        "aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000_blood_normal.read1.fastq.gz",
        "aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000_blood_normal.read2.fastq.gz",
        "aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000_blood_normal.vcf",
        "aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000_blood_tumor.read1.fastq.gz",
        "aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000_blood_tumor.read2.fastq.gz",
        "aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000_blood_tumor.vcf",
        "bbbbbbbb11111111bbbbbbbb11111111bbbbbbbb11111111bbbbbbbb11111111_blood_normal.read1.fastq.gz",
        "bbbbbbbb11111111bbbbbbbb11111111bbbbbbbb11111111bbbbbbbb11111111_blood_normal.read2.fastq.gz",
        "bbbbbbbb11111111bbbbbbbb11111111bbbbbbbb11111111bbbbbbbb11111111_blood_normal.vcf",
        "target_regions.bed",
    ]:
        expected_checksum = calculate_sha256(submission_dir / "files" / file)
        observed_checksum = calculate_sha256(working_dir_path / "files" / file)

        assert expected_checksum == observed_checksum


def are_dir_trees_equal(dir1, dir2):
    """
    Compare two directories recursively. Files in each directory are
    assumed to be equal if their names and contents are equal.

    :param dir1: First directory path
    :param dir2: Second directory path
    :return: True if the directory trees are the same and
        there were no errors while accessing the directories or files,
        False otherwise.
    """
    dirs_cmp = filecmp.dircmp(dir1, dir2)
    if len(dirs_cmp.left_only) > 0 or len(dirs_cmp.right_only) > 0 or len(dirs_cmp.funny_files) > 0:
        return False
    (_, mismatch, errors) = filecmp.cmpfiles(dir1, dir2, dirs_cmp.common_files, shallow=False)
    if len(mismatch) > 0 or len(errors) > 0:
        return False
    for common_dir in dirs_cmp.common_dirs:
        new_dir1 = os.path.join(dir1, common_dir)
        new_dir2 = os.path.join(dir2, common_dir)
        if not are_dir_trees_equal(new_dir1, new_dir2):
            return False
    return True


@mock_aws
def test_upload_download_submission(
    working_dir_path,
    tmpdir_factory,
    temp_config_file_path,
    config_model,
    remote_bucket,
):
    submission_dir = Path("tests/mock_files/submissions/valid_submission")
    transaction_id = "aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000"

    shutil.copytree(
        submission_dir / "encrypted_files",
        working_dir_path / "encrypted_files",
        dirs_exist_ok=True,
    )
    shutil.copytree(submission_dir / "metadata", working_dir_path / "metadata", dirs_exist_ok=True)

    # upload encrypted submission
    testargs = [
        "upload",
        "--submission-dir",
        str(working_dir_path),
        "--config-file",
        temp_config_file_path,
    ]
    runner = CliRunner()
    with mock.patch(
        "grz_cli.models.config.S3Options.__getattr__",
        lambda self, name: None if name == "endpoint_url" else AttributeError,
    ):
        result = runner.invoke(grz_cli.cli.cli, testargs, catch_exceptions=False)

    assert result.exit_code == 0, result.output

    # download
    download_dir = tmpdir_factory.mktemp("submission_download")
    download_dir_path = Path(download_dir.strpath)

    # upload encrypted submission
    testargs = [
        "download",
        "--submission-id",
        transaction_id,
        "--output-dir",
        str(download_dir_path),
        "--config-file",
        temp_config_file_path,
    ]
    runner = CliRunner()

    result = runner.invoke(grz_cli.cli.cli, testargs, catch_exceptions=False)

    assert result.exit_code == 0, result.output

    assert are_dir_trees_equal(
        working_dir_path / "encrypted_files",
        download_dir_path / "encrypted_files",
    ), "Encrypted files are different!"
    assert are_dir_trees_equal(
        working_dir_path / "metadata",
        download_dir_path / "metadata",
    ), "Metadata is different!"
