"""Integration tests for the ``grz-cli submit`` command."""

from pathlib import Path

import click.testing
import pytest
import yaml
from grz_cli.cli import build_cli
from pytest_mock import MockerFixture


@pytest.fixture
def stub_subcommands(mocker: MockerFixture):
    """Replace the subcommands ``submit`` invokes with no-op mocks.

    ``ctx.invoke`` accepts any callable when the target is not a click
    ``Command``, so swapping the module-level imports lets the real submit
    callback run end-to-end without executing real validation, encryption,
    or upload.

    :returns: A parent mock recording every subcommand call, so their relative order is
        observable. ``mock_calls`` names them, for example ``call.validate(...)``.
    """
    mocker.patch("grz_cli.commands.submit.check_version_and_exit_if_needed")
    recorder = mocker.MagicMock()
    for name in ("validate", "encrypt", "upload"):
        recorder.attach_mock(mocker.patch(f"grz_cli.commands.submit.{name}"), name)
    return recorder


def _called_subcommands(recorder) -> list[str]:
    """Names of the subcommands ``recorder`` saw, in call order."""
    return [call[0] for call in recorder.mock_calls]


def _subcommand_kwargs(recorder) -> dict[str, dict]:
    """The keyword arguments each subcommand was called with, keyed by subcommand name."""
    return {call[0]: call[2] for call in recorder.mock_calls}


def _invoke_submit(*extra_args: str) -> click.testing.Result:
    runner = click.testing.CliRunner()
    return runner.invoke(build_cli(), ["submit", *extra_args], catch_exceptions=False)


def test_submit_runs_subcommands_in_order(
    s3_config_path: Path,
    submission_dir: Path,
    stub_subcommands,
) -> None:
    """Encrypting before validating, or uploading before encrypting, would ship bad data."""
    result = _invoke_submit(
        "--submission-dir",
        str(submission_dir),
        "--config-file",
        str(s3_config_path),
    )

    assert result.exit_code == 0, result.output
    assert _called_subcommands(stub_subcommands) == ["validate", "encrypt", "upload"]


def test_submit_hands_each_subcommand_the_submission_it_was_given(
    s3_config_path: Path,
    submission_dir: Path,
    stub_subcommands,
) -> None:
    """All three must act on the same directory, and encrypt must refuse to skip the validation logs."""
    result = _invoke_submit(
        "--submission-dir",
        str(submission_dir),
        "--config-file",
        str(s3_config_path),
    )
    assert result.exit_code == 0, result.output

    kwargs = _subcommand_kwargs(stub_subcommands)
    assert {name: call["submission_dir"] for name, call in kwargs.items()} == {
        "validate": submission_dir,
        "encrypt": submission_dir,
        "upload": submission_dir,
    }
    assert kwargs["encrypt"]["check_validation_logs"] is True, "submit must not encrypt unvalidated data"


def test_submit_accepts_multiple_config_files(
    tmp_path: Path,
    s3_config_path: Path,
    submission_dir: Path,
    stub_subcommands,
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
