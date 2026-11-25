"""Tests of grzctl delete check, purge-records and mark-done."""

import datetime
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import boto3
import click.testing
import grzctl.cli
import pytest
from grz_db.models.submission import ChangeRequestEnum, SubmissionDb, SubmissionStateEnum
from grz_pydantic_models.common import as_aware_datetime
from moto import mock_aws

SUBMITTER_ID = "000000000"
SUBMISSION_ID = f"{SUBMITTER_ID}_2025-01-01_00000000"
OTHER_ID = f"{SUBMITTER_ID}_2025-01-01_11111111"
TAN_G = "3f9a" + "0" * 56 + "0c21"
OTHER_TAN_G = "8b1e" + "1" * 56 + "00d4"


@pytest.fixture
def s3() -> Iterator[Any]:
    """The test config's archives and inbox in a mocked S3, empty."""
    with mock_aws():
        client = boto3.client("s3", region_name="us-east-1")
        for bucket in ("consented", "non_consented", "inbox"):
            client.create_bucket(Bucket=bucket)
        yield client


def _keys(s3: Any, bucket: str) -> set[str]:
    return {o["Key"] for o in s3.list_objects_v2(Bucket=bucket).get("Contents", [])}


def _archive(s3: Any, submission_id: str) -> None:
    for key in ("files/a.fastq.gz.c4gh", "metadata/metadata.json", "logs/progress_upload.cjson"):
        s3.put_object(Bucket="consented", Key=f"{submission_id}/{key}", Body=b"x" * 10)


def _add(db: SubmissionDb, submission_id: str = SUBMISSION_ID, tan_g: str = TAN_G) -> str:
    db.add_submission(submission_id)
    db.modify_submission(submission_id, "tan_g", tan_g)
    db.modify_submission(submission_id, "local_case_id", "patient-1")
    return submission_id


def _request_change(db: SubmissionDb, submission_id: str, instructed_at: str = "2023-01-05") -> None:
    db.add_change_request(
        submission_id,
        ChangeRequestEnum.DELETE,
        requester_name="Dr. Erika Mustermann",
        requester_email="erika.mustermann@example.org",
        requested_at=datetime.date(2023, 1, 4),
        request_email_content="Bitte löschen.",
        data={"instructed_at": instructed_at, "submitter_id": SUBMITTER_ID},
    )


def _deleted(db: SubmissionDb, submission_id: str, on: datetime.date | None = None) -> None:
    at = datetime.datetime.combine(on, datetime.time(), datetime.UTC) if on else None
    db.update_submission_state(submission_id, SubmissionStateEnum.DELETING, timestamp=at)
    db.purge_personal_data(submission_id)
    db.update_submission_state(submission_id, SubmissionStateEnum.DELETED, timestamp=at)


def _states(db: SubmissionDb, submission_id: str) -> list[SubmissionStateEnum]:
    return [s.state for s in sorted(db.get_submission(submission_id).states, key=lambda s: (s.timestamp, s.id))]


def _flat(text: str) -> str:
    """*text* with its line breaks and runs of spaces, which a table's wrapping inserts, made single spaces."""
    return " ".join(text.split())


def _grzctl(config_path: Path, *args: str, input: str | None = None) -> click.testing.Result:
    return click.testing.CliRunner().invoke(
        grzctl.cli.build_cli(), ["--config", str(config_path), "delete", *args], input=input
    )


class TestCheck:
    def test_reports_nothing_when_every_deletion_is_complete(
        self, s3: Any, db: SubmissionDb, migrated_database_config_path: Path
    ):
        _add(db)
        _deleted(db, SUBMISSION_ID)

        result = _grzctl(migrated_database_config_path, "check")

        assert result.exit_code == 0, result.output
        assert "deleted entirely" in result.stderr

    def test_finds_and_repairs_restored_storage(self, s3: Any, db: SubmissionDb, migrated_database_config_path: Path):
        _add(db)
        _deleted(db, SUBMISSION_ID)
        deleted_at = db.get_submission(SUBMISSION_ID).deleted_at
        _archive(s3, SUBMISSION_ID)

        found = _grzctl(migrated_database_config_path, "check")

        assert found.exit_code != 0
        assert "storage restored" in found.stdout
        assert "the consented archive holds 2 objects again" in _flat(found.stdout)

        repaired = _grzctl(migrated_database_config_path, "check", "--repair")

        assert repaired.exit_code == 0, repaired.output
        assert _keys(s3, "consented") == {f"{SUBMISSION_ID}/logs/progress_upload.cjson", f"{SUBMISSION_ID}/deleted"}
        assert _states(db, SUBMISSION_ID)[-1] == SubmissionStateEnum.DELETED
        assert db.get_submission(SUBMISSION_ID).deleted_at == deleted_at
        assert db.get_submission(SUBMISSION_ID).get_latest_state().data["repeated"] == "storage restored from a backup"
        assert _grzctl(migrated_database_config_path, "check").exit_code == 0

    def test_finds_and_repairs_a_restored_database(
        self, s3: Any, db: SubmissionDb, migrated_database_config_path: Path
    ):
        """A tombstone in an archive, but the database, restored to before the deletion, records none."""
        _add(db)
        _archive(s3, SUBMISSION_ID)
        s3.put_object(Bucket="consented", Key=f"{SUBMISSION_ID}/deleted", Body=b"")

        found = _grzctl(migrated_database_config_path, "check")

        assert found.exit_code != 0
        assert "database restored" in found.stdout

        repaired = _grzctl(migrated_database_config_path, "check", "--repair")

        assert repaired.exit_code == 0, repaired.output
        submission = db.get_submission(SUBMISSION_ID)
        assert submission.is_deleted
        assert submission.local_case_id == f"deleted-{SUBMISSION_ID}"
        assert "database restored" in submission.get_latest_state().data["repeated"]
        assert _keys(s3, "consented") == {f"{SUBMISSION_ID}/logs/progress_upload.cjson", f"{SUBMISSION_ID}/deleted"}

    def test_reports_but_does_not_repair_an_interrupted_deletion_or_an_unknown_tombstone(
        self, s3: Any, db: SubmissionDb, migrated_database_config_path: Path
    ):
        _add(db)
        db.update_submission_state(SUBMISSION_ID, SubmissionStateEnum.DELETING)
        s3.put_object(Bucket="consented", Key=f"{OTHER_ID}/deleted", Body=b"")

        result = _grzctl(migrated_database_config_path, "check", "--repair")

        assert result.exit_code != 0
        assert "deletion interrupted" in result.stdout
        assert "not in the database" in result.stdout
        assert _states(db, SUBMISSION_ID) == [SubmissionStateEnum.DELETING]


class TestPurgeRecords:
    @pytest.fixture
    def deletions(self, s3: Any, db: SubmissionDb) -> None:
        """One deletion three years before 2026-01-10, one a day later."""
        for submission_id, tan_g, deleted_on in (
            (SUBMISSION_ID, TAN_G, datetime.date(2023, 1, 10)),
            (OTHER_ID, OTHER_TAN_G, datetime.date(2023, 1, 11)),
        ):
            _add(db, submission_id, tan_g)
            _request_change(db, submission_id)
            _deleted(db, submission_id, on=deleted_on)

    def test_dry_run_lists_what_is_due(self, deletions, db: SubmissionDb, migrated_database_config_path: Path):
        result = _grzctl(migrated_database_config_path, "purge-records", "--due-before", "2026-01-10", "--dry-run")

        assert result.exit_code == 0, result.output
        assert SUBMISSION_ID in result.stdout and OTHER_ID not in result.stdout
        assert "2023-01-10" in result.stdout and "2026-01-10" in result.stdout and "2023-01-05" in result.stdout
        assert db.get_submission(SUBMISSION_ID).tan_g == TAN_G

    def test_purges_what_is_due_and_names_the_mailbox_folders(
        self, deletions, db: SubmissionDb, migrated_database_config_path: Path
    ):
        result = _grzctl(
            migrated_database_config_path, "purge-records", "--due-before", "2026-01-10", "--yes-i-really-mean-it"
        )

        assert result.exit_code == 0, result.output
        assert db.get_submission(SUBMISSION_ID).tan_g is None
        assert db.get_submission(OTHER_ID).tan_g == OTHER_TAN_G
        assert [s.id for s in db.list_change_requests()] == [OTHER_ID]
        assert "2023-01-05  000000000" in result.stderr

    def test_nothing_due(self, deletions, migrated_database_config_path: Path):
        result = _grzctl(migrated_database_config_path, "purge-records", "--due-before", "2026-01-09")

        assert result.exit_code == 0, result.output
        assert "No records are due by 2026-01-09" in result.stderr


class TestMarkDone:
    def test_records_a_deletion_by_hand_on_its_date(
        self, s3: Any, db: SubmissionDb, migrated_database_config_path: Path
    ):
        _add(db)
        _request_change(db, SUBMISSION_ID)
        s3.put_object(Bucket="consented", Key=f"{SUBMISSION_ID}/logs/progress_upload.cjson", Body=b"log")
        s3.put_object(Bucket="inbox", Key=f"{SUBMISSION_ID}/metadata/metadata.json", Body=b"")
        s3.put_object(Bucket="inbox", Key=f"{SUBMISSION_ID}/cleaned", Body=b"")

        result = _grzctl(migrated_database_config_path, "mark-done", SUBMISSION_ID, "--deleted-on", "2026-03-02")

        assert result.exit_code == 0, result.output
        submission = db.get_submission(SUBMISSION_ID)
        assert submission.is_deleted
        assert as_aware_datetime(submission.deleted_at) == datetime.datetime(2026, 3, 2, tzinfo=datetime.UTC)
        assert submission.local_case_id == f"deleted-{SUBMISSION_ID}"
        assert submission.get_latest_state().data["deleted_by_hand"] is True

    @pytest.mark.parametrize(
        ("setup", "message"),
        [
            pytest.param("no-change-request", "has no Delete change request", id="no-change-request"),
            pytest.param("data-left", "still holds 2 objects", id="data-left"),
            pytest.param("deleted", "is recorded already", id="deleted"),
        ],
    )
    def test_refuses(self, s3: Any, db: SubmissionDb, migrated_database_config_path: Path, setup: str, message: str):
        _add(db)
        if setup != "no-change-request":
            _request_change(db, SUBMISSION_ID)
        if setup == "data-left":
            _archive(s3, SUBMISSION_ID)
        if setup == "deleted":
            _deleted(db, SUBMISSION_ID)

        result = _grzctl(migrated_database_config_path, "mark-done", SUBMISSION_ID, "--deleted-on", "2026-03-02")

        assert result.exit_code != 0
        assert message in result.output

    def test_refuses_a_date_in_the_future(self, s3: Any, db: SubmissionDb, migrated_database_config_path: Path):
        _add(db)
        tomorrow = (datetime.date.today() + datetime.timedelta(days=1)).isoformat()

        result = _grzctl(migrated_database_config_path, "mark-done", SUBMISSION_ID, "--deleted-on", tomorrow)

        assert result.exit_code != 0
        assert "in the future" in result.output
