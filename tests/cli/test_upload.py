import filecmp
import os
import shutil
from importlib.metadata import version
from pathlib import Path
from unittest import mock

import grz_cli.cli
import grzctl.cli
import pytest
from click.testing import CliRunner
from grz_common.progress import EncryptionState, FileProgressLogger
from grz_common.workers.submission import Submission, SubmissionValidationError


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


def test_upload_download_submission(
    working_dir_path,
    tmpdir_factory,
    remote_bucket_with_version,
    temp_s3_db_config_file_path,
    initiated_db_test_connection,  # necessary to initiate DB
):
    submission_dir = Path("tests/mock_files/submissions/valid_submission")
    env = {"GRZ_DB__AUTHOR__PRIVATE_KEY_PASSPHRASE": "test"}

    shutil.copytree(submission_dir / "files", working_dir_path / "files", dirs_exist_ok=True)
    shutil.copytree(
        submission_dir / "encrypted_files",
        working_dir_path / "encrypted_files",
        dirs_exist_ok=True,
    )
    shutil.copytree(submission_dir / "metadata", working_dir_path / "metadata", dirs_exist_ok=True)

    logs_dir = working_dir_path / "logs"
    logs_dir.mkdir()
    progress_file = logs_dir / "progress_encrypt.cjson"
    submission = Submission(
        metadata_dir=working_dir_path / "metadata",
        files_dir=working_dir_path / "files",
    )
    progress_logger = FileProgressLogger[EncryptionState](progress_file)
    for file_path, file_metadata in submission.files.items():
        progress_logger.set_state(
            file_path,
            file_metadata,
            state=EncryptionState(encryption_successful=True),
        )

    with mock.patch(
        "grz_common.models.s3.S3Options.__getattr__",
        lambda self, name: None if name == "endpoint_url" else AttributeError,
    ):
        # upload encrypted submission
        upload_args = [
            "upload",
            "--submission-dir",
            str(working_dir_path),
            "--config-file",
            temp_s3_db_config_file_path,
        ]

        runner = CliRunner(env=env)
        cli = grz_cli.cli.build_cli()
        result = runner.invoke(cli, upload_args, catch_exceptions=False)

        assert result.exit_code == 0, result.output
        assert len(result.output) != 0, result.stderr

        submission_id = result.stdout.strip()

        objects_in_bucket = {obj.key: obj for obj in remote_bucket_with_version.objects.all()}
        assert len(objects_in_bucket) > 0, "Upload failed: No objects were found in the mock S3 bucket!"

        assert objects_in_bucket[f"{submission_id}/version"].get()["Body"].read().decode("utf-8") == version("grz-cli")

        cli = grzctl.cli.build_cli()

        add_args = [
            "db",
            "--config-file",
            str(temp_s3_db_config_file_path),
            "submission",
            "add",
            submission_id,
        ]
        runner.invoke(cli, add_args, catch_exceptions=False)

        # download
        download_dir = tmpdir_factory.mktemp("submission_download")
        download_dir_path = Path(download_dir.strpath)

        # download encrypted submission
        download_args = [
            "download",
            "--submission-id",
            submission_id,
            "--output-dir",
            str(download_dir_path),
            "--config-file",
            str(temp_s3_db_config_file_path),
            "--no-update-db",
            "--populate",
        ]
        result = runner.invoke(cli, download_args, catch_exceptions=False)

        assert result.exit_code == 0, result.output

    assert are_dir_trees_equal(
        working_dir_path / "encrypted_files",
        download_dir_path / "encrypted_files",
    ), "Encrypted files are different!"
    assert are_dir_trees_equal(
        working_dir_path / "metadata",
        download_dir_path / "metadata",
    ), "Metadata is different!"


def test_upload_aborts_on_incomplete_encryption(working_dir_path, temp_s3_config_file_path, remote_bucket_with_version):
    """Verify that the upload command fails if the encryption log marks a file as not successful."""
    submission_dir = Path("tests/mock_files/submissions/valid_submission")
    shutil.copytree(submission_dir / "files", working_dir_path / "files", dirs_exist_ok=True)
    shutil.copytree(submission_dir / "metadata", working_dir_path / "metadata", dirs_exist_ok=True)
    (working_dir_path / "encrypted_files").mkdir()

    logs_dir = working_dir_path / "logs"
    logs_dir.mkdir()
    progress_file = logs_dir / "progress_encrypt.cjson"
    submission = Submission(
        metadata_dir=working_dir_path / "metadata",
        files_dir=working_dir_path / "files",
    )
    progress_logger = FileProgressLogger[EncryptionState](progress_file)

    files_iter = iter(submission.files.items())

    # Mark the first file as having failed encryption
    failed_file_path, failed_file_metadata = next(files_iter)
    progress_logger.set_state(
        failed_file_path,
        failed_file_metadata,
        state=EncryptionState(encryption_successful=False, errors=["Interrupted"]),
    )

    # Mark the rest as successful
    for file_path, file_metadata in files_iter:
        progress_logger.set_state(
            file_path,
            file_metadata,
            state=EncryptionState(encryption_successful=True),
        )

    # Attempt upload
    upload_args = [
        "upload",
        "--submission-dir",
        str(working_dir_path),
        "--config-file",
        temp_s3_config_file_path,
    ]
    runner = CliRunner()
    cli = grz_cli.cli.build_cli()
    result = runner.invoke(cli, upload_args, catch_exceptions=True)

    assert result.exit_code != 0
    assert isinstance(result.exc_info[1], SubmissionValidationError)
    error_message = str(result.exc_info[1])
    assert "Will not upload" in error_message
    relative_failed_path = failed_file_path.relative_to(working_dir_path / "files")
    assert str(relative_failed_path) in error_message

    # Ensure it really did fail
    objects_in_bucket = list(remote_bucket_with_version.objects.all())
    assert len(objects_in_bucket) == 1, "Upload should not have happened!"


def test_upload_aborts_if_encryption_log_missing(
    working_dir_path, temp_s3_config_file_path, remote_bucket_with_version
):
    """Verify that the upload command fails if the encryption log is missing entirely."""
    submission_dir = Path("tests/mock_files/submissions/valid_submission")
    shutil.copytree(submission_dir / "files", working_dir_path / "files", dirs_exist_ok=True)
    shutil.copytree(submission_dir / "metadata", working_dir_path / "metadata", dirs_exist_ok=True)
    (working_dir_path / "encrypted_files").mkdir()
    (working_dir_path / "logs").mkdir()

    # Attempt upload
    upload_args = [
        "upload",
        "--submission-dir",
        str(working_dir_path),
        "--config-file",
        temp_s3_config_file_path,
    ]
    runner = CliRunner()
    cli = grz_cli.cli.build_cli()
    result = runner.invoke(cli, upload_args, catch_exceptions=True)

    assert result.exit_code != 0
    assert isinstance(result.exc_info[1], SubmissionValidationError)
    error_message = str(result.exc_info[1])
    assert "Will not upload" in error_message

    # Check if at least one of the files is listed as unencrypted
    assert "target_regions.bed" in error_message

    # Ensure it really did fail
    objects_in_bucket = list(remote_bucket_with_version.objects.all())
    assert len(objects_in_bucket) == 1, "Upload should not have happened!"


def test_upload_workflow_succeeds_with_symlink_in_files_dir(
    working_dir_path,
    temp_identifiers_config_file_path,
    temp_keys_config_file_path,
    temp_s3_config_file_path,
    remote_bucket_with_version,
):
    """
    End-to-end smoke test for LE submissions where `files/` contains symlinks.

    We intentionally use `--no-grz-check` here to exercise the Python fallback
    validation path in a deterministic way.
    """
    submission_dir = Path("tests/mock_files/submissions/valid_submission")

    shutil.copytree(submission_dir / "files", working_dir_path / "files", dirs_exist_ok=True)
    shutil.copytree(submission_dir / "metadata", working_dir_path / "metadata", dirs_exist_ok=True)

    # Replace one submission file with a symlink pointing to data outside the submission directory.
    symlink_path = working_dir_path / "files" / "target_regions.bed"
    external_dir = working_dir_path / "real_data"
    external_dir.mkdir(parents=True, exist_ok=True)
    external_target = external_dir / symlink_path.name

    shutil.move(symlink_path, external_target)
    try:
        symlink_path.symlink_to(external_target)
    except (OSError, NotImplementedError) as e:
        pytest.skip(f"Symlinks not supported in this environment: {e}")

    runner = CliRunner()
    cli = grz_cli.cli.build_cli()
    config_args = [
        "--config-file",
        str(temp_identifiers_config_file_path),
        "--config-file",
        str(temp_keys_config_file_path),
        "--config-file",
        str(temp_s3_config_file_path),
    ]

    # validate (fallback)
    result = runner.invoke(
        cli,
        [
            "validate",
            "--submission-dir",
            str(working_dir_path),
            "--no-grz-check",
            *config_args,
        ],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.output

    # encrypt (must read the symlinked input file)
    result = runner.invoke(
        cli,
        [
            "encrypt",
            "--submission-dir",
            str(working_dir_path),
            *config_args,
        ],
        catch_exceptions=False,
    )
    assert result.exit_code == 0, result.output

    # upload (must not fail due to symlink usage earlier)
    with mock.patch(
        "grz_common.models.s3.S3Options.__getattr__",
        lambda self, name: None if name == "endpoint_url" else AttributeError,
    ):
        result = runner.invoke(
            cli,
            [
                "upload",
                "--submission-dir",
                str(working_dir_path),
                *config_args,
            ],
            catch_exceptions=False,
        )
    assert result.exit_code == 0, result.output

    submission_id = result.stdout.strip()
    assert submission_id, "Expected `grz-cli upload` to output a submission id."

    objects_in_bucket = {obj.key for obj in remote_bucket_with_version.objects.all()}
    assert f"{submission_id}/metadata/metadata.json" in objects_in_bucket
    assert f"{submission_id}/files/target_regions.bed.c4gh" in objects_in_bucket
