import errno
import re

import grz_cli.cli
import grzctl.cli
import pytest
import yaml
from click.testing import CliRunner
from grz_common.exceptions import ConfigurationError, EncryptionError, IncompleteSubmissionError
from grz_common.progress import FileProgressLogger, ValidationState
from grz_common.utils.checksums import calculate_sha256
from grz_common.utils.crypt import Crypt4GH
from grz_common.workers.submission import Submission

from .common import SUBMISSION_DIR, copy_submission


def test_encrypt_submission_protect_overwrite(
    working_dir_path,
    temp_keys_config_file_path,
    tmpdir_factory: pytest.TempdirFactory,
):
    copy_submission(working_dir_path, "files", "metadata")

    testargs = [
        "encrypt",
        "--submission-dir",
        str(working_dir_path),
        "--config-file",
        temp_keys_config_file_path,
        "--no-check-validation-logs",
    ]

    runner = CliRunner()
    cli = grz_cli.cli.build_cli()
    # run encrypt once to build logs
    runner.invoke(cli, testargs, catch_exceptions=False)

    # running again should not error because cache is used
    runner.invoke(cli, testargs, catch_exceptions=False)

    # removing the cache and running again should error without force
    (working_dir_path / "logs" / "progress_encrypt.cjson").unlink()
    with pytest.raises(EncryptionError, match=re.escape("already exists. Delete it or use --force to overwrite it.")):
        runner.invoke(cli, testargs, catch_exceptions=False)


def test_encrypt_with_an_unreadable_public_key_fails_as_a_configuration_error(
    working_dir_path, keys_config_content, tmp_path
):
    """A key that the setup cannot provide is not the file's encryption failing."""
    copy_submission(working_dir_path, "files", "metadata")
    not_a_key = tmp_path / "not_a_key.pub"
    not_a_key.write_text("not a key")
    keys_config_content["keys"]["grz_public_key_path"] = str(not_a_key)
    config_file = tmp_path / "config.keys.yaml"
    config_file.write_text(yaml.dump(keys_config_content))

    testargs = [
        "encrypt",
        "--submission-dir",
        str(working_dir_path),
        "--config-file",
        str(config_file),
        "--no-check-validation-logs",
    ]
    result = CliRunner().invoke(grz_cli.cli.build_cli(), testargs)

    assert isinstance(result.exception, ConfigurationError), result.output


def test_decrypt_submission(working_dir_path, temp_grzctl_keys_config_file_path):
    copy_submission(working_dir_path, "encrypted_files", "metadata")

    testargs = [
        "--config",
        temp_grzctl_keys_config_file_path,
        "decrypt",
        "--submission-dir",
        str(working_dir_path),
        "--no-update-db",
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
        expected_checksum = calculate_sha256(SUBMISSION_DIR / "files" / file)
        observed_checksum = calculate_sha256(working_dir_path / "files" / file)

        assert expected_checksum == observed_checksum


def test_decrypt_lets_an_error_the_file_did_not_cause_through(
    working_dir_path, temp_grzctl_keys_config_file_path, monkeypatch
):
    """A full disk is no fault of the submitter, so decryption does not turn it into a decryption error."""
    copy_submission(working_dir_path, "encrypted_files", "metadata")

    def fail(*_args, **_kwargs):
        raise OSError(errno.ENOSPC, "No space left on device")

    monkeypatch.setattr(Crypt4GH, "decrypt_file", fail)
    testargs = [
        "--config",
        temp_grzctl_keys_config_file_path,
        "decrypt",
        "--submission-dir",
        str(working_dir_path),
        "--no-update-db",
    ]

    result = CliRunner().invoke(grzctl.cli.build_cli(), testargs)

    assert isinstance(result.exception, OSError), result.output


def test_encrypt_decrypt_submission(
    working_dir_path,
    temp_keys_config_file_path,
    temp_grzctl_keys_config_file_path,
    # crypt4gh_grz_private_key_file_path,
    tmpdir_factory: pytest.TempdirFactory,
):
    copy_submission(working_dir_path, "files", "metadata")

    # first, encrypt the data
    testargs = [
        "encrypt",
        "--submission-dir",
        str(working_dir_path),
        "--config-file",
        temp_keys_config_file_path,
        "--no-check-validation-logs",
    ]

    runner = CliRunner()
    cli = grz_cli.cli.build_cli()
    result = runner.invoke(cli, testargs, catch_exceptions=False)

    assert result.exit_code == 0, result.output

    # then, decrypt the data again
    testargs = [
        "--config",
        temp_grzctl_keys_config_file_path,
        "decrypt",
        "--submission-dir",
        str(working_dir_path),
        "--no-update-db",
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
        expected_checksum = calculate_sha256(SUBMISSION_DIR / "files" / file)
        observed_checksum = calculate_sha256(working_dir_path / "files" / file)

        assert expected_checksum == observed_checksum


def test_encrypt_succeeds_with_valid_logs(working_dir_path, temp_keys_config_file_path):
    """Verify that the encrypt command succeeds if validation logs are present and mark all files as valid."""
    copy_submission(working_dir_path, "files", "metadata")
    logs_dir = working_dir_path / "logs"
    logs_dir.mkdir()

    # Create valid validation logs
    submission = Submission(
        metadata_dir=working_dir_path / "metadata",
        files_dir=working_dir_path / "files",
    )
    checksum_progress_logger = FileProgressLogger[ValidationState](logs_dir / "progress_checksum_validation.cjson")
    seq_data_progress_logger = FileProgressLogger[ValidationState](
        logs_dir / "progress_sequencing_data_validation.cjson"
    )

    for file_path, file_metadata in submission.files.items():
        checksum_progress_logger.set_state(
            file_path,
            file_metadata,
            state=ValidationState(validation_passed=True),
        )
        if file_metadata.file_type in ("fastq", "bam"):
            seq_data_progress_logger.set_state(
                file_path,
                file_metadata,
                state=ValidationState(validation_passed=True),
            )

    # Attempt encryption
    encrypt_args = [
        "encrypt",
        "--submission-dir",
        str(working_dir_path),
        "--config-file",
        temp_keys_config_file_path,
    ]
    runner = CliRunner()
    cli = grz_cli.cli.build_cli()
    result = runner.invoke(cli, encrypt_args, catch_exceptions=False)

    assert result.exit_code == 0, result.output
    assert (working_dir_path / "encrypted_files").exists()
    # Check if at least one file got encrypted
    assert len(list((working_dir_path / "encrypted_files").iterdir())) > 0


def test_encrypt_aborts_on_incomplete_validation(working_dir_path, temp_keys_config_file_path):
    """Verify that the encrypt command fails if the validation log marks a file as not successful."""
    copy_submission(working_dir_path, "files", "metadata")
    logs_dir = working_dir_path / "logs"
    logs_dir.mkdir()

    # Create validation logs with one failed file
    submission = Submission(
        metadata_dir=working_dir_path / "metadata",
        files_dir=working_dir_path / "files",
    )
    checksum_progress_logger = FileProgressLogger[ValidationState](logs_dir / "progress_checksum_validation.cjson")
    files_iter = iter(submission.files.items())
    failed_file_path, failed_file_metadata = next(files_iter)

    checksum_progress_logger.set_state(
        failed_file_path,
        failed_file_metadata,
        state=ValidationState(validation_passed=False, errors=["Checksum mismatch"]),
    )

    # Mark the rest as successful
    for file_path, file_metadata in files_iter:
        checksum_progress_logger.set_state(
            file_path,
            file_metadata,
            state=ValidationState(validation_passed=True),
        )

    # Attempt encryption
    encrypt_args = [
        "encrypt",
        "--submission-dir",
        str(working_dir_path),
        "--config-file",
        temp_keys_config_file_path,
    ]
    runner = CliRunner()
    cli = grz_cli.cli.build_cli()
    result = runner.invoke(cli, encrypt_args, catch_exceptions=True)

    assert result.exit_code != 0
    assert isinstance(result.exc_info[1], IncompleteSubmissionError)
    error_message = str(result.exc_info[1])
    assert "Will not encrypt" in error_message
    assert str(failed_file_path) in error_message
    assert not (working_dir_path / "encrypted_files").exists()


def test_encrypt_aborts_if_validation_log_missing(working_dir_path, temp_keys_config_file_path):
    """Verify that the encrypt command fails if the validation log is missing entirely."""
    copy_submission(working_dir_path, "files", "metadata")
    (working_dir_path / "logs").mkdir()  # create empty logs dir

    # Attempt encryption
    encrypt_args = [
        "encrypt",
        "--submission-dir",
        str(working_dir_path),
        "--config-file",
        temp_keys_config_file_path,
    ]
    runner = CliRunner()
    cli = grz_cli.cli.build_cli()
    result = runner.invoke(cli, encrypt_args, catch_exceptions=True)

    assert result.exit_code != 0
    assert isinstance(result.exc_info[1], IncompleteSubmissionError)
    error_message = str(result.exc_info[1])
    assert "Will not encrypt" in error_message

    # Check if at least one of the files is listed as unvalidated
    assert "target_regions.bed" in error_message
    assert not (working_dir_path / "encrypted_files").exists()
