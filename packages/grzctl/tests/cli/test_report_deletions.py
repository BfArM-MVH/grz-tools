"""How the reports count deletions, and what they warn about them."""

import csv
import datetime
import logging
from pathlib import Path

import pytest
from grz_db.models.submission import ChangeRequestEnum, SubmissionDb, SubmissionStateEnum
from grzctl.commands.report import (
    _dump_overview_report,
    _next_purge_date,
    _warn_about_open_deletions,
    _warn_about_records_due,
)
from grzctl.models.config import GrzctlConfig

SUBMITTER_ID = "260914050"
GRZ_ID = "GRZK00007"


def _add(db: SubmissionDb, suffix: str, uploaded: datetime.date, tan_g: str | None = None) -> str:
    submission_id = f"{SUBMITTER_ID}_{uploaded.isoformat()}_{suffix}"
    db.add_submission(submission_id)
    db.modify_submission(submission_id, "submitter_id", SUBMITTER_ID)
    db.modify_submission(submission_id, "data_node_id", GRZ_ID)
    db.modify_submission(submission_id, "submission_uploaded_date", uploaded)
    if tan_g is not None:
        db.modify_submission(submission_id, "tan_g", tan_g)
    return submission_id


def _request_deletion(db: SubmissionDb, submission_id: str) -> None:
    db.add_change_request(
        submission_id,
        ChangeRequestEnum.DELETE,
        requester_name="Dr. Erika Mustermann",
        requester_email="erika.mustermann@example.org",
        requested_at=datetime.date(2026, 1, 5),
        request_email_content="Bitte löschen.",
    )


def _delete(db: SubmissionDb, submission_id: str, at: datetime.datetime) -> None:
    for state in (SubmissionStateEnum.DELETING, SubmissionStateEnum.DELETED):
        db.update_submission_state(submission_id, state, timestamp=at)


def _deletions(db: SubmissionDb, tmp_path: Path, year: int, quarter: int) -> int:
    path = tmp_path / f"overview_{year}_{quarter}.tsv"
    _dump_overview_report(path, db, year, quarter)
    with open(path, newline="", encoding="utf-8") as overview:
        rows = [row for row in csv.DictReader(overview, delimiter="\t") if row["submitterId"] == SUBMITTER_ID]
    return int(rows[0]["number_of_deletions"]) if rows else 0


class TestNumberOfDeletions:
    def test_counts_a_deletion_in_the_quarter_it_was_carried_out(self, db: SubmissionDb, tmp_path: Path):
        """Uploaded in Q1, deleted at the very end of Q2: listed in the report for Q2, as the BfArM SOP asks."""
        submission_id = _add(db, "00000000", datetime.date(2026, 2, 1))
        _request_deletion(db, submission_id)
        _delete(db, submission_id, datetime.datetime(2026, 6, 30, 23, 59, 59, tzinfo=datetime.UTC))

        assert [_deletions(db, tmp_path, 2026, quarter) for quarter in (1, 2, 3)] == [0, 1, 0]

    def test_a_request_alone_is_no_deletion(self, db: SubmissionDb, tmp_path: Path):
        submission_id = _add(db, "00000000", datetime.date(2026, 2, 1))
        _request_deletion(db, submission_id)

        assert _deletions(db, tmp_path, 2026, 1) == 0

    def test_an_interrupted_deletion_is_no_deletion(self, db: SubmissionDb, tmp_path: Path):
        submission_id = _add(db, "00000000", datetime.date(2026, 2, 1))
        db.update_submission_state(
            submission_id, SubmissionStateEnum.DELETING, timestamp=datetime.datetime(2026, 2, 3, tzinfo=datetime.UTC)
        )

        assert _deletions(db, tmp_path, 2026, 1) == 0

    def test_a_repeated_deletion_counts_once_in_the_quarter_of_the_first(self, db: SubmissionDb, tmp_path: Path):
        """Grzctl delete check --repair records Deleted again; the deletion stays in its first quarter."""
        submission_id = _add(db, "00000000", datetime.date(2026, 2, 1))
        _delete(db, submission_id, datetime.datetime(2026, 2, 3, tzinfo=datetime.UTC))
        db.update_submission_state(
            submission_id, SubmissionStateEnum.DELETED, timestamp=datetime.datetime(2026, 5, 3, tzinfo=datetime.UTC)
        )

        assert [_deletions(db, tmp_path, 2026, quarter) for quarter in (1, 2)] == [1, 0]


class TestWarnings:
    def test_warns_about_a_requested_deletion_never_carried_out(
        self, db: SubmissionDb, caplog: pytest.LogCaptureFixture
    ):
        open_request = _add(db, "00000000", datetime.date(2026, 2, 1))
        _request_deletion(db, open_request)
        carried_out = _add(db, "11111111", datetime.date(2026, 2, 1))
        _request_deletion(db, carried_out)
        _delete(db, carried_out, datetime.datetime(2026, 2, 3, tzinfo=datetime.UTC))

        with caplog.at_level(logging.WARNING):
            _warn_about_open_deletions(db)

        assert "1 submission(s) have a Delete change request but were never deleted" in caplog.text
        assert open_request in caplog.text and carried_out not in caplog.text

    def test_no_warning_without_open_deletions(self, db: SubmissionDb, caplog: pytest.LogCaptureFixture):
        _add(db, "00000000", datetime.date(2026, 2, 1))

        with caplog.at_level(logging.WARNING):
            _warn_about_open_deletions(db)

        assert caplog.text == ""

    @pytest.mark.parametrize(
        ("today", "due"),
        [
            pytest.param(datetime.date(2027, 1, 15), datetime.date(2027, 6, 30), id="january"),
            pytest.param(datetime.date(2027, 4, 15), datetime.date(2027, 9, 30), id="april"),
            pytest.param(datetime.date(2027, 10, 15), datetime.date(2028, 3, 31), id="october"),
            pytest.param(datetime.date(2027, 12, 31), datetime.date(2028, 3, 31), id="new-years-eve"),
        ],
    )
    def test_next_purge_date_is_the_end_of_the_quarter_after_the_current_one(
        self, today: datetime.date, due: datetime.date
    ):
        assert _next_purge_date(today) == due

    def test_warns_about_records_due_before_the_next_report(
        self, db: SubmissionDb, migrated_database_config: GrzctlConfig, caplog: pytest.LogCaptureFixture
    ):
        """Run on 15.01.2027, it covers everything due by 30.06.2027: deleted by 30.06.2024."""
        due = _add(db, "00000000", datetime.date(2024, 1, 2), tan_g="a" * 64)
        _delete(db, due, datetime.datetime(2024, 6, 30, 12, tzinfo=datetime.UTC))
        later = _add(db, "11111111", datetime.date(2024, 1, 2), tan_g="b" * 64)
        _delete(db, later, datetime.datetime(2024, 7, 1, 12, tzinfo=datetime.UTC))

        with caplog.at_level(logging.WARNING):
            _warn_about_records_due(migrated_database_config.db.database_url, datetime.date(2027, 1, 15))

        assert "The records of 1 deletion(s) fall due by 30.06.2027" in caplog.text
        assert "grzctl delete purge-records --due-before 2027-06-30" in caplog.text
