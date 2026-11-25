import datetime
from pathlib import Path

import click
import click.testing
import grzctl.cli
import pytest
import yaml
from grz_db.models.submission import SubmissionDb, SubmissionStateEnum, SubmissionStateLog
from grzctl.commands.delete.confirmation import build_confirmation
from grzctl.models.deletion_request import DeletionRequest

GRZ_ID = "GRZK00007"
SUBMITTER_ID = "000000000"
TAN_G_A = "3f9a" + "0" * 56 + "0c21"
TAN_G_B = "8b1e" + "1" * 56 + "00d4"
TAN_G_UNKNOWN = "c" * 64
TODAY = datetime.date(2026, 12, 1)


def _request(*tan_gs: str, instructed_at: datetime.date = datetime.date(2026, 11, 21)) -> DeletionRequest:
    return DeletionRequest(
        instructed_at=instructed_at,
        submitter_id=SUBMITTER_ID,
        requester_name="Dr. Erika Mustermann",
        requester_email="erika.mustermann@example.org",
        requested_at=datetime.date(2026, 11, 20),
        reason="Widerruf der Einwilligung",
        data_node_id=GRZ_ID,
        tan_g=list(tan_gs),
    )


def _add(db: SubmissionDb, suffix: str, tan_g: str, *states: tuple[SubmissionStateEnum, datetime.datetime]) -> str:
    """A submission with *tan_g* and *states*, each written with its own timestamp."""
    submission_id = f"{SUBMITTER_ID}_2025-01-01_{suffix}"
    db.add_submission(submission_id)
    db.modify_submission(submission_id, "tan_g", tan_g)
    with db.transaction() as session:
        for state, timestamp in states:
            session.add(
                SubmissionStateLog(
                    submission_id=submission_id,
                    state=state,
                    timestamp=timestamp,
                    author_name="alice",
                    signature="unsigned",
                )
            )
    return submission_id


def _deleted_on(day: datetime.date) -> tuple[tuple[SubmissionStateEnum, datetime.datetime], ...]:
    at = datetime.datetime.combine(day, datetime.time(10), datetime.UTC)
    return (SubmissionStateEnum.DELETING, at), (SubmissionStateEnum.DELETED, at)


def test_confirms_each_tan_g_in_full(db: SubmissionDb):
    _add(db, "00000000", TAN_G_A, *_deleted_on(datetime.date(2026, 11, 25)))

    text = build_confirmation(_request(TAN_G_A, TAN_G_UNKNOWN), db, GRZ_ID, TODAY)

    assert text.startswith("Betreff: Löschung von Datensätzen im genomDE MV\n")
    assert "hiermit bestätigen wir die Löschung der Datensätze gemäß Ihrer Anweisung vom 21.11.2026" in text
    assert "Leistungserbringers 000000000 vom 20.11.2026" in text
    assert f"  {TAN_G_A}  gelöscht am 25.11.2026\n" in text
    assert f"  {TAN_G_UNKNOWN}  liegt beim GRZ GRZK00007 nicht vor\n" in text
    assert "spätestens bis zum 23.02.2027 gelöscht" in text
    assert "im Quartalsbericht für Q4/2026 in der Spalte number_of_deletions" in text
    assert "bis längstens 25.11.2029 aufbewahrt" in text


def test_names_each_quarter_of_a_deletion_across_a_quarter_boundary(db: SubmissionDb):
    _add(db, "00000000", TAN_G_A, *_deleted_on(datetime.date(2026, 12, 31)))
    _add(db, "11111111", TAN_G_B, *_deleted_on(datetime.date(2027, 1, 2)))

    text = build_confirmation(_request(TAN_G_A, TAN_G_B), db, GRZ_ID, datetime.date(2027, 1, 5))

    assert "Quartalsbericht für Q4/2026 und Q1/2027" in text
    assert "spätestens bis zum 02.04.2027 gelöscht" in text
    assert "bis längstens 02.01.2030 aufbewahrt" in text


def test_only_tan_gs_not_present(db: SubmissionDb):
    text = build_confirmation(_request(TAN_G_UNKNOWN), db, GRZ_ID, TODAY)

    assert "zu Ihrer Anweisung vom 21.11.2026" in text
    assert "bestätigen" not in text
    assert "Backups" not in text and "Quartalsbericht" not in text


@pytest.mark.parametrize(
    "states",
    [
        pytest.param((), id="not-deleted"),
        pytest.param(
            ((SubmissionStateEnum.DELETING, datetime.datetime(2026, 11, 25, tzinfo=datetime.UTC)),), id="interrupted"
        ),
    ],
)
def test_refuses_before_every_submission_is_deleted(db: SubmissionDb, states):
    _add(db, "00000000", TAN_G_A, *_deleted_on(datetime.date(2026, 11, 25)))
    unfinished = _add(db, "11111111", TAN_G_B, *states)

    with pytest.raises(click.ClickException, match=f"Not all datasets are deleted yet: {unfinished}"):
        build_confirmation(_request(TAN_G_A, TAN_G_B), db, GRZ_ID, TODAY)


def test_refuses_a_request_whose_records_are_purged(db: SubmissionDb):
    with pytest.raises(click.ClickException, match="Its records are purged"):
        build_confirmation(_request(TAN_G_A, instructed_at=datetime.date(2023, 11, 21)), db, GRZ_ID, TODAY)


def test_command_prints_the_confirmation(db: SubmissionDb, migrated_database_config_path: Path, tmp_path: Path):
    _add(db, "00000000", TAN_G_A, *_deleted_on(datetime.date(2026, 11, 25)))
    data_file = tmp_path / "antrag.yaml"
    data_file.write_text(yaml.safe_dump(_request(TAN_G_A).model_dump(mode="json"), allow_unicode=True))

    result = click.testing.CliRunner().invoke(
        grzctl.cli.build_cli(),
        ["--config", str(migrated_database_config_path), "delete", "confirmation", "--data-file", str(data_file)],
    )

    assert result.exit_code == 0, result.output
    assert f"  {TAN_G_A}  gelöscht am 25.11.2026\n" in result.stdout
