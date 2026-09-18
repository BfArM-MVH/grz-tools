import signal

import click.testing
import grzctl.cli
import pytest


def test_help():
    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()
    result = runner.invoke(cli, ["--help"])
    assert result.exit_code == 0, result.stderr


def test_sigterm_stops_a_run_like_ctrl_c():
    """A KeyboardInterrupt is what makes the running step record ``interrupted``."""
    previous = signal.getsignal(signal.SIGTERM)
    try:
        grzctl.cli._stop_on_sigterm()
        with pytest.raises(KeyboardInterrupt):
            signal.raise_signal(signal.SIGTERM)
    finally:
        signal.signal(signal.SIGTERM, previous)
