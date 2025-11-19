import filecmp
import json
import os
import shutil
from importlib.metadata import version
from pathlib import Path
from unittest import mock

import grz_cli.cli
import grzctl.cli
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
    temp_s3_config_file_path,
    remote_bucket,
):
    submission_dir = Path("tests/mock_files/submissions/valid_submission")

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
            temp_s3_config_file_path,
        ]

        runner = CliRunner()
        cli = grz_cli.cli.build_cli()
        result = runner.invoke(cli, upload_args, catch_exceptions=False)

        assert result.exit_code == 0, result.output
        assert len(result.output) != 0, result.stderr

        submission_id = result.stdout.strip()

        objects_in_bucket = {obj.key: obj for obj in remote_bucket.objects.all()}
        assert len(objects_in_bucket) > 0, "Upload failed: No objects were found in the mock S3 bucket!"

        assert objects_in_bucket[f"{submission_id}/version"].get()["Body"].read().decode("utf-8") == version("grz-cli")

        # Get the upload date of the metadata file from S3
        metadata_s3_object = objects_in_bucket[f"{submission_id}/metadata/metadata.json"]
        upload_date = metadata_s3_object.last_modified.date()

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
            temp_s3_config_file_path,
        ]
        cli = grzctl.cli.build_cli()
        result = runner.invoke(cli, download_args, catch_exceptions=False)

        assert result.exit_code == 0, result.output

    assert are_dir_trees_equal(
        working_dir_path / "encrypted_files",
        download_dir_path / "encrypted_files",
    ), "Encrypted files are different!"

    # The submission date in the downloaded metadata is updated, so we expect the metadata to be different.
    with open(download_dir_path / "metadata" / "metadata.json") as f:
        downloaded_metadata = json.load(f)

    assert downloaded_metadata["submission"]["submissionDate"] == upload_date.isoformat()


def test_upload_aborts_on_incomplete_encryption(
    working_dir_path,
    temp_s3_config_file_path,
    remote_bucket,
):
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
    objects_in_bucket = list(remote_bucket.objects.all())
    assert len(objects_in_bucket) == 0, "Upload should not have happened!"


def test_upload_aborts_if_encryption_log_missing(
    working_dir_path,
    temp_s3_config_file_path,
    remote_bucket,
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
    objects_in_bucket = list(remote_bucket.objects.all())
    assert len(objects_in_bucket) == 0, "Upload should not have happened!"
