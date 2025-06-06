"""
Tests for the Prüfbericht submission functionality.
"""

import importlib.resources
import shutil

import click.testing
import grzctl

from .. import mock_files


def test_archive(temp_s3_config_file_path, remote_bucket, working_dir_path):
    submission_dir_ptr = importlib.resources.files(mock_files).joinpath("submissions", "valid_submission")
    with importlib.resources.as_file(submission_dir_ptr) as submission_dir:
        shutil.copytree(submission_dir / "encrypted_files", working_dir_path / "encrypted_files", dirs_exist_ok=True)
        shutil.copytree(submission_dir / "metadata", working_dir_path / "metadata", dirs_exist_ok=True)

        args = [
            "archive",
            "--config-file",
            temp_s3_config_file_path,
            "--submission-dir",
            str(working_dir_path),
        ]

        runner = click.testing.CliRunner(
            env={
                "GRZ_S3_OPTIONS__ENDPOINT_URL": "",
            }
        )
        cli = grzctl.cli.build_cli()
        result = runner.invoke(cli, args, catch_exceptions=False)

    assert result.exit_code == 0, result.output
