"""How commands outside grzctl delete treat a submission whose deletion has started."""

import logging
from pathlib import Path

import click
import click.testing
import grzctl.cli
import pytest
from grz_db.errors import SubmissionDeletedError
from grz_db.models.submission import SubmissionDb, SubmissionStateEnum
from grzctl.dbcontext import DbContext
from grzctl.models.config import GrzctlConfig

SUBMISSION_ID = "000000000_2025-01-01_00000000"


@pytest.mark.parametrize("state", [SubmissionStateEnum.DELETING, SubmissionStateEnum.DELETED])
def test_db_context_passes_the_refusal_on(db: SubmissionDb, migrated_database_config: GrzctlConfig, state):
    """Not as 'Failed to connect to DB', which would send one looking for a connection problem."""
    db.add_submission(SUBMISSION_ID)
    db.update_submission_state(SUBMISSION_ID, SubmissionStateEnum.REPORTED)
    db.update_submission_state(SUBMISSION_ID, state)

    with (
        pytest.raises(SubmissionDeletedError),
        DbContext(
            configuration=migrated_database_config,
            submission_id=SUBMISSION_ID,
            start_state=SubmissionStateEnum.CLEANING,
            end_state=SubmissionStateEnum.CLEANED,
        ),
    ):
        pass


def test_main_reports_the_refusal_without_a_traceback(
    monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
):
    @click.command()
    def refused():
        raise SubmissionDeletedError(SUBMISSION_ID)

    monkeypatch.setattr(grzctl.cli, "build_cli", lambda: refused)
    monkeypatch.setattr("sys.argv", ["grzctl"])

    with caplog.at_level(logging.ERROR), pytest.raises(SystemExit) as exit_info:
        grzctl.cli.main()

    assert exit_info.value.code == 1
    assert f"Submission '{SUBMISSION_ID}' is being or has been deleted" in caplog.text
    assert "Traceback" not in caplog.text


@pytest.mark.parametrize("state", ["Deleting", "deleted"])
def test_db_submission_update_refuses_the_deletion_states(
    db: SubmissionDb, migrated_database_config_path: Path, state: str
):
    """Set by hand, they would record a deletion that never happened, which the signed log cannot take back."""
    db.add_submission(SUBMISSION_ID)

    result = click.testing.CliRunner().invoke(
        grzctl.cli.build_cli(),
        ["--config", str(migrated_database_config_path), "db", "submission", "update", SUBMISSION_ID, state],
    )

    assert result.exit_code != 0
    assert "grzctl delete run" in result.stderr
    assert db.get_submission(SUBMISSION_ID).states == []
