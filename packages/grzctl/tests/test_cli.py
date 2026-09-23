import signal

import click
import click.testing
import grzctl.cli
import pytest
from grz_common.exceptions import ConfigurationError


def test_help():
    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()
    result = runner.invoke(cli, ["--help"])
    assert result.exit_code == 0, result.stderr


def test_main_reports_a_grz_error_without_traceback(monkeypatch, caplog):
    """An expected failure ends with its message and exit code 1."""

    @click.command()
    def failing():
        raise ConfigurationError("the key cannot be read")

    monkeypatch.setattr(grzctl.cli, "build_cli", lambda: failing)
    monkeypatch.setattr("sys.argv", ["grzctl"])

    with pytest.raises(SystemExit) as exit_info:
        grzctl.cli.main()

    assert exit_info.value.code == 1
    assert "the key cannot be read" in caplog.text


def test_sigterm_stops_a_run_like_ctrl_c():
    """A KeyboardInterrupt is what makes the running step record ``interrupted``."""
    previous = signal.getsignal(signal.SIGTERM)
    try:
        grzctl.cli._stop_on_sigterm()
        with pytest.raises(KeyboardInterrupt):
            signal.raise_signal(signal.SIGTERM)
    finally:
        signal.signal(signal.SIGTERM, previous)
