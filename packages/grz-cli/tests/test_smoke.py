def test_import_cli():
    import grz_cli.cli  # noqa: F401


def test_help(capsys):
    import click.testing

    from grz_cli.cli import build_cli

    result = click.testing.CliRunner().invoke(build_cli(), ["--help"])
    assert result.exit_code == 0
    assert "Usage:" in result.output
