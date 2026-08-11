"""
Tests for the grzctl list functionality.
"""

import importlib.resources
import json
from unittest import mock

import click.testing
import grzctl
from grz_common.progress import EncryptionState, FileProgressLogger
from grz_common.workers.submission import Submission

from .. import mock_files
from .common import copy_submission


def test_list(temp_s3_config_file_path, remote_bucket_with_version, working_dir_path, tmp_path):
    submission_dir_ptr = importlib.resources.files(mock_files).joinpath("submissions", "valid_submission")
    with importlib.resources.as_file(submission_dir_ptr) as submission_dir:
        copy_submission(working_dir_path, "files", "encrypted_files", "metadata", source=submission_dir)

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


def test_list_with_partial_env(temp_s3_config_file_path, remote_bucket_with_version, working_dir_path, tmp_path):
    """If database configuration is partially-populated via environment variables, it should still be ignored."""
    submission_dir_ptr = importlib.resources.files(mock_files).joinpath("submissions", "valid_submission")
    with importlib.resources.as_file(submission_dir_ptr) as submission_dir:
        copy_submission(working_dir_path, "files", "encrypted_files", "metadata", source=submission_dir)

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
