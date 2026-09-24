import contextlib
import errno
import re
from collections.abc import Iterator
from pathlib import Path
from unittest.mock import patch

import crypt4gh.header
import crypt4gh.keys
import grz_cli.cli
import grz_common.exceptions as grzexc
import grzctl.cli
import pytest
import yaml
from click.testing import CliRunner
from grz_common.progress import FileProgressLogger, ValidationState
from grz_common.utils.checksums import calculate_sha256
from grz_common.utils.crypt import Crypt4GH
from grz_common.workers.submission import Submission

from ..conftest import _grzctl_config_dict
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
    with pytest.raises(
        grzexc.EncryptionError, match=re.escape("already exists. Delete it or use --force to overwrite it.")
    ):
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

    assert isinstance(result.exception, grzexc.ConfigurationError), result.output


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


def test_encrypt_signs_with_the_submitter_key(
    working_dir_path,
    temp_keys_config_file_path,
    crypt4gh_grz_private_key_file_path,
    crypt4gh_submitter_public_key_file_path,
    crypt4gh_grz_public_key_file_path,
):
    """The header of an encrypted file names the submitter's key as the sender, not a random key."""
    copy_submission(working_dir_path, "files", "metadata")
    testargs = [
        "encrypt",
        "--submission-dir",
        str(working_dir_path),
        "--config-file",
        temp_keys_config_file_path,
        "--no-check-validation-logs",
    ]
    result = CliRunner().invoke(grz_cli.cli.build_cli(), testargs, catch_exceptions=False)
    assert result.exit_code == 0, result.output

    # crypt4gh works with the raw 32 bytes of each key
    keys = [(0, Crypt4GH.retrieve_private_key(crypt4gh_grz_private_key_file_path).private_bytes_raw(), None)]
    submitter_public_key = Crypt4GH.retrieve_public_key(crypt4gh_submitter_public_key_file_path).public_bytes_raw()
    other_public_key = Crypt4GH.retrieve_public_key(crypt4gh_grz_public_key_file_path).public_bytes_raw()
    encrypted_file = next((working_dir_path / "encrypted_files").rglob("*.c4gh"))
    with open(encrypted_file, "rb") as f:
        crypt4gh.header.deconstruct(f, keys, sender_pubkey=submitter_public_key)
    with open(encrypted_file, "rb") as f, pytest.raises(ValueError):
        crypt4gh.header.deconstruct(f, keys, sender_pubkey=other_public_key)


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
    assert isinstance(result.exc_info[1], grzexc.IncompleteSubmissionError)
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
    assert isinstance(result.exc_info[1], grzexc.IncompleteSubmissionError)
    error_message = str(result.exc_info[1])
    assert "Will not encrypt" in error_message

    # Check if at least one of the files is listed as unvalidated
    assert "target_regions.bed" in error_message
    assert not (working_dir_path / "encrypted_files").exists()


GRZ_PRIVATE_KEY = str(Path("tests/mock_files/grz_mock_private_key.sec").resolve())
SUBMITTER_PRIVATE_KEY = str(Path("tests/mock_files/submitter_mock_private_key.sec").resolve())
SUBMITTER_PUBLIC_KEY = Path("tests/mock_files/submitter_mock_public_key.pub")
CONSENTED_PUBLIC_KEY = str(Path("tests/mock_files/archive_consented.pub").resolve())
CONSENTED_PRIVATE_KEY = Path("tests/mock_files/archive_consented.sec").resolve()
NON_CONSENTED_PUBLIC_KEY = str(Path("tests/mock_files/archive_non_consented.pub").resolve())
NON_CONSENTED_PRIVATE_KEY = str(Path("tests/mock_files/archive_non_consented.sec").resolve())
SUBMITTER_ID = "260914050"
"""The submitter named in the example submission's metadata."""
DECRYPTED_FILE = "target_regions.bed"


def _write_grzctl_config(tmp_path: Path, leistungserbringer: dict) -> Path:
    """Write a grzctl config with the given inboxes, and with a public key of its own for each archive."""
    config = _grzctl_config_dict(leistungserbringer=leistungserbringer)
    config["archives"]["consented"]["public_key_path"] = CONSENTED_PUBLIC_KEY
    config["archives"]["non_consented"]["public_key_path"] = NON_CONSENTED_PUBLIC_KEY
    config_path = tmp_path / "config.grzctl.yaml"
    config_path.write_text(yaml.safe_dump(config))
    return config_path


def _run_grzctl(config_path: Path, command: str, submission_dir: Path, *args: str):
    runner = CliRunner()
    cli = grzctl.cli.build_cli()
    return runner.invoke(
        cli, ["--config", str(config_path), command, "--submission-dir", str(submission_dir), "--no-update-db", *args]
    )


@contextlib.contextmanager
def _loaded_key_paths() -> Iterator[list[str]]:
    """Record the path of every private key file that is loaded."""
    loaded: list[str] = []
    retrieve_private_key = Crypt4GH.retrieve_private_key

    def _record(seckey_path, *args, **kwargs):
        loaded.append(str(seckey_path))
        return retrieve_private_key(seckey_path, *args, **kwargs)

    with patch.object(Crypt4GH, "retrieve_private_key", side_effect=_record):
        yield loaded


def _assert_decrypted(working_dir_path: Path):
    expected_checksum = calculate_sha256(SUBMISSION_DIR / "files" / DECRYPTED_FILE)
    assert calculate_sha256(working_dir_path / "files" / DECRYPTED_FILE) == expected_checksum


def test_decrypt_uses_the_inbox_key_of_the_submitter_in_the_metadata(working_dir_path, tmp_path):
    """Only the inboxes of the metadata's submitter count, and the key that they share loads once."""
    copy_submission(working_dir_path, "encrypted_files", "metadata")
    config_path = _write_grzctl_config(
        tmp_path,
        {
            "000000000": {"inbox_buckets": {"inbox": {"private_key_path": SUBMITTER_PRIVATE_KEY}}},
            SUBMITTER_ID: {
                "inbox_buckets": {
                    "first": {"private_key_path": GRZ_PRIVATE_KEY},
                    "second": {"private_key_path": GRZ_PRIVATE_KEY},
                }
            },
        },
    )

    with _loaded_key_paths() as loaded:
        result = _run_grzctl(config_path, "decrypt", working_dir_path)

    assert result.exit_code == 0, result.output
    assert loaded == [GRZ_PRIVATE_KEY]
    _assert_decrypted(working_dir_path)


def test_grzctl_encrypt_encrypts_for_the_archive_and_signs_with_the_signing_key(working_dir_path, tmp_path):
    """``grzctl encrypt`` re-encrypts the files for the archive that the consent selects, here the consented one.

    grzctl has no archive private key, so the test decrypts a file with the archive's key itself.
    """
    copy_submission(working_dir_path, "files", "metadata")
    config_path = _write_grzctl_config(
        tmp_path, {SUBMITTER_ID: {"inbox_buckets": {"testing": {"private_key_path": GRZ_PRIVATE_KEY}}}}
    )
    config = yaml.safe_load(config_path.read_text())
    config["archives"]["signing_key_path"] = SUBMITTER_PRIVATE_KEY
    config_path.write_text(yaml.safe_dump(config))

    result = _run_grzctl(config_path, "encrypt", working_dir_path, "--no-check-validation-logs")
    assert result.exit_code == 0, result.output

    encrypted_file_path = working_dir_path / "encrypted_files" / f"{DECRYPTED_FILE}.c4gh"
    with open(encrypted_file_path, "rb") as encrypted_file:
        packet = next(crypt4gh.header.parse(encrypted_file))
    # an X25519 header packet starts with the 4 bytes of the method, then the sender's public key
    assert packet[4:36] == crypt4gh.keys.get_public_key(SUBMITTER_PUBLIC_KEY), "the signing key signs the files"

    decrypted_file_path = tmp_path / DECRYPTED_FILE
    Crypt4GH.decrypt_file(
        encrypted_file_path, decrypted_file_path, Crypt4GH.retrieve_private_key(CONSENTED_PRIVATE_KEY)
    )
    assert calculate_sha256(decrypted_file_path) == calculate_sha256(SUBMISSION_DIR / "files" / DECRYPTED_FILE)


def test_decrypt_fails_if_the_inbox_key_does_not_open_the_files(working_dir_path, tmp_path):
    """The error says that the key does not fit, but never shows the key material."""
    copy_submission(working_dir_path, "encrypted_files", "metadata")
    config_path = _write_grzctl_config(
        tmp_path, {SUBMITTER_ID: {"inbox_buckets": {"testing": {"private_key": CONSENTED_PRIVATE_KEY.read_text()}}}}
    )

    result = _run_grzctl(config_path, "decrypt", working_dir_path)

    assert isinstance(result.exception, grzexc.DecryptionError), result.output
    message = str(result.exception)
    assert "the private key does not open its Crypt4GH header" in message
    for line in CONSENTED_PRIVATE_KEY.read_text().splitlines():
        assert line not in message


def test_decrypt_fails_for_a_submitter_missing_from_the_config(working_dir_path, tmp_path):
    """Without an inbox of the submitter, no key is loaded, so another submitter's key cannot stand in."""
    copy_submission(working_dir_path, "encrypted_files", "metadata")
    config_path = _write_grzctl_config(
        tmp_path, {"000000000": {"inbox_buckets": {"inbox": {"private_key_path": GRZ_PRIVATE_KEY}}}}
    )

    with _loaded_key_paths() as loaded:
        result = _run_grzctl(config_path, "decrypt", working_dir_path)

    assert isinstance(result.exception, grzexc.ConfigurationError), result.output
    assert f"Submitter '{SUBMITTER_ID}' is not in the config" in str(result.exception)
    assert loaded == []


def test_decrypt_fails_if_the_inboxes_of_the_submitter_use_different_keys(working_dir_path, tmp_path):
    copy_submission(working_dir_path, "encrypted_files", "metadata")
    config_path = _write_grzctl_config(
        tmp_path,
        {
            SUBMITTER_ID: {
                "inbox_buckets": {
                    "a": {"private_key_path": GRZ_PRIVATE_KEY},
                    "b": {"private_key_path": NON_CONSENTED_PRIVATE_KEY},
                }
            }
        },
    )

    with _loaded_key_paths() as loaded:
        result = _run_grzctl(config_path, "decrypt", working_dir_path)

    assert isinstance(result.exception, grzexc.ConfigurationError), result.output
    message = str(result.exception)
    assert f"Key 1: leistungserbringer.{SUBMITTER_ID}.inbox_buckets.a.private_key_path." in message
    assert f"Key 2: leistungserbringer.{SUBMITTER_ID}.inbox_buckets.b.private_key_path." in message
    assert loaded == []
    assert not (working_dir_path / "files" / DECRYPTED_FILE).exists()
