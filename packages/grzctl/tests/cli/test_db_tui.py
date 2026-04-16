import logging

import click.testing
import grzctl.cli


def test_db_tui_passes_quarter_year(monkeypatch, blank_database_config_path):
    captured: dict[str, object] = {}

    class DummyBrowser:
        def __init__(self, *, database, public_keys, quarter=None, year=None, **kwargs):
            captured["quarter"] = quarter
            captured["year"] = year

        def run(self):
            captured["ran"] = True

    import grzctl.commands.db.cli as db_cli

    monkeypatch.setattr(db_cli, "DatabaseBrowser", DummyBrowser)

    root_logger = logging.getLogger()
    old_handlers = list(root_logger.handlers)
    try:
        runner = click.testing.CliRunner(catch_exceptions=False)
        cli = grzctl.cli.build_cli()
        result = runner.invoke(
            cli,
            [
                "db",
                "--config-file",
                str(blank_database_config_path),
                "tui",
                "--quarter",
                "2",
                "--year",
                "2025",
            ],
        )
        assert result.exit_code == 0, result.stderr
        assert captured["quarter"] == 2
        assert captured["year"] == 2025
        assert captured.get("ran") is True
    finally:
        root_logger.handlers[:] = old_handlers


def test_db_tui_defaults_to_current_when_omitted(monkeypatch, blank_database_config_path):
    captured: dict[str, object] = {}

    class DummyBrowser:
        def __init__(self, *, database, public_keys, quarter=None, year=None, **kwargs):
            captured["quarter"] = quarter
            captured["year"] = year

        def run(self):
            captured["ran"] = True

    import grzctl.commands.db.cli as db_cli

    monkeypatch.setattr(db_cli, "DatabaseBrowser", DummyBrowser)

    root_logger = logging.getLogger()
    old_handlers = list(root_logger.handlers)
    try:
        runner = click.testing.CliRunner(catch_exceptions=False)
        cli = grzctl.cli.build_cli()
        result = runner.invoke(
            cli,
            [
                "db",
                "--config-file",
                str(blank_database_config_path),
                "tui",
            ],
        )
        assert result.exit_code == 0, result.stderr
        assert captured["quarter"] is None
        assert captured["year"] is None
        assert captured.get("ran") is True
    finally:
        root_logger.handlers[:] = old_handlers
