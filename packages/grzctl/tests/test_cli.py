import click.testing
import grzctl.cli


def test_help():
    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()
    result = runner.invoke(cli, ["--help"])
    assert result.exit_code == 0, result.stderr
