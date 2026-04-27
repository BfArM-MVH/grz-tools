import datetime
import logging

import click.testing
import grzctl.cli
import grzctl.commands.db.tui as db_tui
import pytest


@pytest.fixture
def captured_browser_args(monkeypatch):
    captured: dict[str, object] = {}

    class DummyBrowser:
        def __init__(self, *, database, public_keys, quarter=None, year=None, **kwargs):
            captured["quarter"] = quarter
            captured["year"] = year

        def run(self):
            captured["ran"] = True

    import grzctl.commands.db.cli as db_cli

    monkeypatch.setattr(db_cli, "DatabaseBrowser", DummyBrowser)
    return captured


def test_db_tui_passes_quarter_year(captured_browser_args, blank_database_config_path):
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
        assert captured_browser_args["quarter"] == 2
        assert captured_browser_args["year"] == 2025
        assert captured_browser_args.get("ran") is True
    finally:
        root_logger.handlers[:] = old_handlers


def test_db_tui_passes_none_when_quarter_year_omitted(captured_browser_args, blank_database_config_path):
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
        assert captured_browser_args["quarter"] is None
        assert captured_browser_args["year"] is None
        assert captured_browser_args.get("ran") is True
    finally:
        root_logger.handlers[:] = old_handlers


def test_resolve_quarter_year_defaults_to_current_when_omitted(monkeypatch):
    def fake_date_to_quarter_year(_today):
        return (3, 2026)

    monkeypatch.setattr(db_tui, "date_to_quarter_year", fake_date_to_quarter_year)

    resolved_quarter, resolved_year = db_tui._resolve_quarter_year(
        datetime.date(2026, 8, 12),
        quarter=None,
        year=None,
    )

    assert resolved_quarter == 3
    assert resolved_year == 2026


def test_resolve_quarter_year_prefers_explicit_values(monkeypatch):
    def fake_date_to_quarter_year(_today):
        return (3, 2026)

    monkeypatch.setattr(db_tui, "date_to_quarter_year", fake_date_to_quarter_year)

    resolved_quarter, resolved_year = db_tui._resolve_quarter_year(
        datetime.date(2026, 8, 12),
        quarter=2,
        year=2025,
    )

    assert resolved_quarter == 2
    assert resolved_year == 2025
