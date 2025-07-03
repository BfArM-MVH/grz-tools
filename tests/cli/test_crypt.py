import shutil
from pathlib import Path

import grz_cli.cli
import grzctl.cli
import pytest
from click.testing import CliRunner
from grz_common.utils.checksums import calculate_sha256


def test_encrypt_submission(
    working_dir_path,
    temp_keys_config_file_path,
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
        temp_keys_config_file_path,
    ]

    runner = CliRunner()
    cli = grz_cli.cli.build_cli()
    result = runner.invoke(cli, testargs, catch_exceptions=False)

    assert result.exit_code == 0, result.output


def test_encrypt_submission_protect_overwrite(
    working_dir_path,
    temp_keys_config_file_path,
    tmpdir_factory: pytest.TempdirFactory,
):
    submission_dir = Path("tests/mock_files/submissions/valid_submission")

    shutil.copytree(submission_dir / "files", working_dir_path / "files", dirs_exist_ok=True)
    shutil.copytree(submission_dir / "metadata", working_dir_path / "metadata", dirs_exist_ok=True)

    testargs = [
        "encrypt",
        "--submission-dir",
        str(working_dir_path),
        "--config-file",
        temp_keys_config_file_path,
    ]

    runner = CliRunner()
    cli = grz_cli.cli.build_cli()
    # run encrypt once to build logs
    runner.invoke(cli, testargs, catch_exceptions=False)

    # running again should not error because cache is used
    runner.invoke(cli, testargs, catch_exceptions=False)

    # removing the cache and running again should error without force
    (working_dir_path / "logs" / "progress_encrypt.cjson").unlink()
    with pytest.raises(RuntimeError, match="already exists. Delete it or use --force to overwrite it."):
        runner.invoke(cli, testargs, catch_exceptions=False)


def test_decrypt_submission(working_dir_path, temp_keys_config_file_path):
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
        temp_keys_config_file_path,
    ]
    runner = CliRunner()
    cli = grzctl.cli.build_cli()
    result = runner.invoke(cli, testargs, catch_exceptions=False)

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


def test_encrypt_decrypt_submission(
    working_dir_path,
    temp_keys_config_file_path,
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
        temp_keys_config_file_path,
    ]

    runner = CliRunner()
    cli = grz_cli.cli.build_cli()
    result = runner.invoke(cli, testargs, catch_exceptions=False)

    assert result.exit_code == 0, result.output

    # then, decrypt the data again
    testargs = [
        "decrypt",
        "--submission-dir",
        str(working_dir_path),
        "--config-file",
        temp_keys_config_file_path,
    ]

    runner = CliRunner()
    cli = grzctl.cli.build_cli()
    result = runner.invoke(cli, testargs, catch_exceptions=False)

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
