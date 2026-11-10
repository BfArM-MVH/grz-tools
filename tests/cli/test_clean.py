"""
Tests for the Prüfbericht submission functionality.
"""

import importlib.resources
import json

import botocore.client
import click.testing
import grz_cli.cli
import grz_common.exceptions as grzexc
import grzctl.cli
import pytest
from botocore.exceptions import ClientError
from grz_common.progress import EncryptionState, FileProgressLogger
from grz_common.workers.submission import Submission
from grzctl.commands.clean import _clean_submission_from_bucket

from .. import mock_files
from .common import copy_submission


def test_clean_and_list(
    temp_s3_config_file_path,
    temp_grzctl_s3_config_file_path,
    temp_grzctl_s3_db_config_file_path,
    remote_bucket_with_version,
    working_dir_path,
    tmp_path,
):
    submission_dir_ptr = importlib.resources.files(mock_files).joinpath("submissions", "valid_submission")
    with importlib.resources.as_file(submission_dir_ptr) as submission_dir:
        copy_submission(working_dir_path, "files", "encrypted_files", "metadata", source=submission_dir)

        # manually set successful encrypted state in progress logs since upload checks for this
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

    # upload encrypted submission
    upload_args = [
        "upload",
        "--submission-dir",
        str(working_dir_path),
        "--config-file",
        temp_s3_config_file_path,
    ]

    runner = click.testing.CliRunner()
    cli = grz_cli.cli.build_cli()
    result_upload = runner.invoke(cli, upload_args, catch_exceptions=False)

    assert result_upload.exit_code == 0, result_upload.output
    assert len(result_upload.output) != 0, result_upload.stderr

    submission_id = result_upload.stdout.strip()

    cli = grzctl.cli.build_cli()

    clean_args = [
        "--config",
        temp_grzctl_s3_config_file_path,
        "clean",
        "--submission-id",
        submission_id,
        "--yes-i-really-mean-it",
        "--no-update-db",
        "--inbox",
        "testing",
    ]

    result_clean = runner.invoke(cli, clean_args, catch_exceptions=False)

    assert result_clean.exit_code == 0, result_clean.output

    uploaded_keys = {o.key for o in remote_bucket_with_version.objects.all()}
    assert len(uploaded_keys) == 3
    assert f"{submission_id}/metadata/metadata.json" in uploaded_keys
    assert f"{submission_id}/cleaned" in uploaded_keys
    assert f"{submission_id}/cleaning" not in uploaded_keys
    # ensure metadata is empty
    assert remote_bucket_with_version.Object(f"{submission_id}/metadata/metadata.json").content_length == 0

    list_args = [
        "--config",
        temp_grzctl_s3_db_config_file_path,
        "list",
        "--json",
        "--show-cleaned",
        "--inbox",
        "testing",
        "--submitter-id",
        "260914050",
    ]

    result_list = runner.invoke(cli, list_args, catch_exceptions=False)

    assert result_list.exit_code == 0, result_list.output

    listed_submissions = json.loads(result_list.stdout.strip())
    assert len(listed_submissions) == 1
    assert listed_submissions[0]["state"] == "cleaned"


def _fail_s3_operation(monkeypatch, operation: str, code: str):
    """Answer every S3 call of *operation* with the error *code*."""
    original_call = botocore.client.BaseClient._make_api_call

    def fail(self, operation_name, kwargs):
        if operation_name == operation:
            raise ClientError({"Error": {"Code": code, "Message": code}}, operation_name)
        return original_call(self, operation_name, kwargs)

    monkeypatch.setattr(botocore.client.BaseClient, "_make_api_call", fail)


def test_clean_reports_rejected_credentials_as_a_configuration_error(s3_config_model, remote_bucket, monkeypatch):
    """Rejected credentials during clean's S3 calls must classify as a configuration error, not a raw ClientError."""
    _fail_s3_operation(monkeypatch, "PutObject", "InvalidAccessKeyId")

    with pytest.raises(grzexc.ConfigurationError):
        _clean_submission_from_bucket(
            s3_config_model.s3.bucket, s3_config_model.s3, "123_2025-01-01_00000000", "'testing'"
        )
