"""
Tests for the grz-cli get-id subcommand.
"""

import importlib.resources

import click.testing
import grz_cli

from .. import mock_files
from .common import copy_submission


def test_get_id(working_dir_path):
    submission_dir_ptr = importlib.resources.files(mock_files).joinpath("submissions", "valid_submission")
    with importlib.resources.as_file(submission_dir_ptr) as submission_dir:
        copy_submission(working_dir_path, "encrypted_files", "metadata", source=submission_dir)

        runner = click.testing.CliRunner()
        cli = grz_cli.cli.build_cli()
        result = runner.invoke(cli, ["get-id", str(submission_dir / "metadata" / "metadata.json")])

        assert result.exit_code == 0, result.output

        submission_id = result.stdout.strip()
        assert submission_id == "260914050_2024-07-15_c64603a7"
