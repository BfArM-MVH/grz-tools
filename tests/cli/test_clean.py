"""
Tests for the Prüfbericht submission functionality.
"""

import importlib.resources
import json
import shutil
from unittest import mock

import click.testing
import grzctl

from .. import mock_files


def test_clean_and_list(temp_s3_config_file_path, remote_bucket, working_dir_path, tmp_path):
    submission_dir_ptr = importlib.resources.files(mock_files).joinpath("submissions", "valid_submission")
    with importlib.resources.as_file(submission_dir_ptr) as submission_dir:
        shutil.copytree(submission_dir / "encrypted_files", working_dir_path / "encrypted_files", dirs_exist_ok=True)
        shutil.copytree(submission_dir / "metadata", working_dir_path / "metadata", dirs_exist_ok=True)

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

        clean_args = [
            "clean",
            "--submission-id",
            submission_id,
            "--config-file",
            temp_s3_config_file_path,
            "--yes-i-really-mean-it",
        ]

        result_clean = runner.invoke(cli, clean_args, catch_exceptions=False)

        assert result_clean.exit_code == 0, result_clean.output

        uploaded_keys = {o.key for o in remote_bucket.objects.all()}
        assert len(uploaded_keys) == 2
        assert f"{submission_id}/metadata/metadata.json" in uploaded_keys
        assert f"{submission_id}/cleaned" in uploaded_keys
        assert f"{submission_id}/cleaning" not in uploaded_keys
        # ensure metadata is empty
        assert remote_bucket.Object(f"{submission_id}/metadata/metadata.json").content_length == 0

        list_args = ["list", "--config-file", temp_s3_config_file_path, "--json", "--show-cleaned"]

        result_list = runner.invoke(cli, list_args, catch_exceptions=False)

        assert result_list.exit_code == 0, result_list.output

        listed_submissions = json.loads(result_list.stdout.strip())
        assert len(listed_submissions) == 1
        assert listed_submissions[0]["state"] == "cleaned"
