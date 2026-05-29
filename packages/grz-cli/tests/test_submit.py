"""Integration tests for the ``grz-cli submit`` command."""

from pathlib import Path

import click.testing
import pytest
import yaml
from grz_cli.cli import build_cli
from pytest_mock import MockerFixture


@pytest.fixture
def stub_subcommands(mocker: MockerFixture) -> dict[str, object]:
    """Replace the subcommands ``submit`` invokes with no-op mocks.

    ``ctx.invoke`` accepts any callable when the target is not a click
    ``Command``, so swapping the module-level imports lets the real submit
    callback run end-to-end without executing real validation, encryption,
    or upload.
    """
    mocker.patch("grz_cli.commands.submit.check_version_and_exit_if_needed")
    return {
        "validate": mocker.patch("grz_cli.commands.submit.validate"),
        "encrypt": mocker.patch("grz_cli.commands.submit.encrypt"),
        "upload": mocker.patch("grz_cli.commands.submit.upload"),
    }


def _invoke_submit(*extra_args: str) -> click.testing.Result:
    runner = click.testing.CliRunner()
    return runner.invoke(build_cli(), ["submit", *extra_args], catch_exceptions=False)


def test_submit_runs_subcommands_in_order(
    s3_config_path: Path,
    submission_dir: Path,
    stub_subcommands: dict[str, object],
) -> None:
    result = _invoke_submit(
        "--submission-dir",
        str(submission_dir),
        "--config-file",
        str(s3_config_path),
    )

    assert result.exit_code == 0, result.output
    stub_subcommands["validate"].assert_called_once()
    stub_subcommands["encrypt"].assert_called_once()
    stub_subcommands["upload"].assert_called_once()


def test_submit_accepts_multiple_config_files(
    tmp_path: Path,
    s3_config_path: Path,
    submission_dir: Path,
    stub_subcommands: dict[str, object],
) -> None:
    """Regression: repeated ``--config-file`` once produced a tuple that crashed config parsing."""
    override_path = tmp_path / "override.yaml"
    override_path.write_text(yaml.dump({"s3": {"access_key": "override-key"}}))

    result = _invoke_submit(
        "--submission-dir",
        str(submission_dir),
        "--config-file",
        str(s3_config_path),
        "--config-file",
        str(override_path),
    )

    assert result.exit_code == 0, result.output
