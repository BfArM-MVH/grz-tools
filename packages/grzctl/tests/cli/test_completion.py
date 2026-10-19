"""grzctl completion prints click's completion script, and db commands can complete submission IDs."""

import socket
import threading
from pathlib import Path

import click.testing
import grzctl.cli
import pytest
from click.shell_completion import ShellComplete

from .conftest import _database_config, _write_config

SUBMISSION_IDS = ["123456789_2024-11-08_d0f805c5", "123456789_2024-11-09_a1b2c3d4", "987654321_2024-11-08_e5f6a7b8"]


@pytest.mark.parametrize("shell", ["bash", "zsh", "fish"])
def test_completion_prints_script(shell: str):
    result = click.testing.CliRunner().invoke(grzctl.cli.build_cli(), ["completion", shell])

    assert result.exit_code == 0, result.stderr
    assert f"_GRZCTL_COMPLETE={shell}_complete" in result.stdout


@pytest.fixture
def config_path_with_submissions(migrated_database_config_path: Path) -> Path:
    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()
    for submission_id in SUBMISSION_IDS:
        result = runner.invoke(
            cli, ["--config", str(migrated_database_config_path), "db", "submission", "add", submission_id]
        )
        assert result.exit_code == 0, result.stderr
    return migrated_database_config_path


def _complete(config_path: Path, args: list[str], incomplete: str) -> list[str]:
    complete = ShellComplete(grzctl.cli.build_cli(), {}, "grzctl", "_GRZCTL_COMPLETE")
    return [item.value for item in complete.get_completions(["--config", str(config_path), *args], incomplete)]


def test_submission_ids_complete_only_when_enabled(config_path_with_submissions: Path, monkeypatch):
    args = ["db", "submission", "show"]
    monkeypatch.delenv("GRZCTL_COMPLETE_SUBMISSION_IDS", raising=False)
    assert _complete(config_path_with_submissions, args, "123") == []

    monkeypatch.setenv("GRZCTL_COMPLETE_SUBMISSION_IDS", "1")
    assert _complete(config_path_with_submissions, args, "123") == SUBMISSION_IDS[:2]
    assert _complete(config_path_with_submissions, ["db", "backfill", "--submission-id"], "987") == SUBMISSION_IDS[2:]


def test_submission_id_completion_swallows_errors(tmp_path: Path, monkeypatch):
    monkeypatch.setenv("GRZCTL_COMPLETE_SUBMISSION_IDS", "1")
    assert _complete(tmp_path / "missing.yaml", ["db", "submission", "show"], "") == []


def test_submission_id_completion_times_out_on_unresponsive_postgresql(tmp_path: Path, monkeypatch):
    """A server that accepts the TCP connection but never answers must not block Tab."""
    monkeypatch.setenv("GRZCTL_COMPLETE_SUBMISSION_IDS", "1")
    with socket.create_server(("127.0.0.1", 0)) as server:
        port = server.getsockname()[1]
        config = _database_config(tmp_path, f"postgresql+psycopg://grzctl@127.0.0.1:{port}/grzctl")
        config_path = _write_config(tmp_path, config)
        completions = []
        # a daemon thread, so a hanging completion fails the test instead of hanging the run
        thread = threading.Thread(
            target=lambda: completions.extend(_complete(config_path, ["db", "submission", "show"], "")),
            daemon=True,
        )
        thread.start()
        thread.join(timeout=10)

        assert not thread.is_alive()
        assert completions == []
