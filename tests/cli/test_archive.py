"""
Tests for the Prüfbericht submission functionality.
"""

import importlib.resources
import shutil

import click.testing

import grz_cli

from .. import mock_files
from .common import working_dir, working_dir_path  # noqa: F401


def test_archive(temp_config_file_path, remote_bucket, working_dir_path):  # noqa: F811
    submission_dir_ptr = importlib.resources.files(mock_files).joinpath("submissions", "valid_submission")
    with importlib.resources.as_file(submission_dir_ptr) as submission_dir:
        shutil.copytree(submission_dir / "encrypted_files", working_dir_path / "encrypted_files", dirs_exist_ok=True)
        shutil.copytree(submission_dir / "metadata", working_dir_path / "metadata", dirs_exist_ok=True)

        args = [
            "archive",
            "--config-file",
            temp_config_file_path,
            "--submission-dir",
            str(working_dir_path),
        ]

        runner = click.testing.CliRunner(
            env={
                "GRZ_S3_OPTIONS__ENDPOINT_URL": "",
            }
        )
        cli = grz_cli.cli.build_cli(grz_mode=True)
        result = runner.invoke(cli, args, catch_exceptions=False)

    assert result.exit_code == 0, result.output
