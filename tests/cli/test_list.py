"""
Tests for the grzctl list functionality.
"""

import importlib.resources
import json
import shutil
from unittest import mock

import click.testing
import grzctl
from grz_common.progress import EncryptionState, FileProgressLogger
from grz_common.workers.submission import Submission

from .. import mock_files


def test_list(temp_s3_config_file_path, remote_bucket, working_dir_path, tmp_path):
    submission_dir_ptr = importlib.resources.files(mock_files).joinpath("submissions", "valid_submission")
    with importlib.resources.as_file(submission_dir_ptr) as submission_dir:
        shutil.copytree(submission_dir / "files", working_dir_path / "files", dirs_exist_ok=True)
        shutil.copytree(submission_dir / "encrypted_files", working_dir_path / "encrypted_files", dirs_exist_ok=True)
        shutil.copytree(submission_dir / "metadata", working_dir_path / "metadata", dirs_exist_ok=True)

        logs_dir = working_dir_path / "logs"
        logs_dir.mkdir(exist_ok=True)
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

        runner = click.testing.CliRunner()
        cli = grzctl.cli.build_cli()
        result_upload = runner.invoke(cli, upload_args, catch_exceptions=False)

        assert result_upload.exit_code == 0, result_upload.output
        assert len(result_upload.output) != 0, result_upload.stderr

        submission_id = result_upload.stdout.strip()

        list_args = ["list", "--config-file", temp_s3_config_file_path, "--json", "--show-cleaned"]

        result_list = runner.invoke(cli, list_args, catch_exceptions=False)

        assert result_list.exit_code == 0, result_list.output

        listed_submissions = json.loads(result_list.stdout.strip())
        assert len(listed_submissions) == 1
        assert listed_submissions[0]["submission_id"] == submission_id
        assert listed_submissions[0]["state"] == "complete"


def test_list_with_partial_env(temp_s3_config_file_path, remote_bucket, working_dir_path, tmp_path):
    """If database configuration is partially-populated via environment variables, it should still be ignored."""
    submission_dir_ptr = importlib.resources.files(mock_files).joinpath("submissions", "valid_submission")
    with importlib.resources.as_file(submission_dir_ptr) as submission_dir:
        shutil.copytree(submission_dir / "files", working_dir_path / "files", dirs_exist_ok=True)
        shutil.copytree(submission_dir / "encrypted_files", working_dir_path / "encrypted_files", dirs_exist_ok=True)
        shutil.copytree(submission_dir / "metadata", working_dir_path / "metadata", dirs_exist_ok=True)

        logs_dir = working_dir_path / "logs"
        logs_dir.mkdir(exist_ok=True)
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

        runner = click.testing.CliRunner()
        cli = grzctl.cli.build_cli()
        result_upload = runner.invoke(cli, upload_args, catch_exceptions=False)

        assert result_upload.exit_code == 0, result_upload.output
        assert len(result_upload.output) != 0, result_upload.stderr

        submission_id = result_upload.stdout.strip()

        list_args = ["list", "--config-file", temp_s3_config_file_path, "--json", "--show-cleaned"]

        result_list = runner.invoke(
            cli, list_args, catch_exceptions=False, env={"GRZ_DB__AUTHOR__PRIVATE_KEY_PASSPHRASE": "secret"}
        )

        assert result_list.exit_code == 0, result_list.output

        listed_submissions = json.loads(result_list.stdout.strip())
        assert len(listed_submissions) == 1
        assert listed_submissions[0]["submission_id"] == submission_id
        assert listed_submissions[0]["state"] == "complete"
