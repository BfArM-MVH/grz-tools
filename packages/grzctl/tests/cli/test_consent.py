"""Tests for ``grzctl consent``."""

import importlib.resources
from pathlib import Path

import click.testing
import grzctl.cli
from grz_pydantic_models_testing.example_metadata import grzctl as grzctl_metadata


def test_consent_details_name_why_a_donor_gives_no_research_consent(tmp_path: Path) -> None:
    """The details table says why a donor gives no research consent."""
    metadata_dir = tmp_path / "metadata"
    metadata_dir.mkdir()
    (metadata_dir / "metadata.json").write_text(
        importlib.resources.files(grzctl_metadata).joinpath("metadata.json").read_text()
    )

    result = click.testing.CliRunner().invoke(
        grzctl.cli.build_cli(),
        ["consent", "--submission-dir", str(tmp_path), "--details", "--date", "2025-09-15"],
        # wide enough that the table does not wrap the reason
        env={"COLUMNS": "300"},
    )

    assert result.exit_code == 0, result.output
    assert "researchConsents[0] has no scope, noScopeJustification 'other patient-related reason'" in result.output
