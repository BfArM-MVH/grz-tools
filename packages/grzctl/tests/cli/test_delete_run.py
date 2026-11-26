import datetime
import json
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import boto3
import click.testing
import grzctl.cli
import pytest
import yaml
from grz_db.models.submission import RequestRawContentType, SubmissionDb, SubmissionStateEnum, SubmissionType
from grzctl.commands.delete import execute
from grzctl.commands.delete.resolve import (
    DeletionOutcome,
    DeletionRequestRefusedError,
    resolve_deletion_request,
)
from grzctl.models.deletion_request import DeletionRequest
from moto import mock_aws

GRZ_ID = "GRZK00007"
SUBMITTER_ID = "000000000"
TAN_G_A = "3f9a" + "0" * 56 + "0c21"
TAN_G_B = "8b1e" + "1" * 56 + "00d4"
TAN_G_UNKNOWN = "c" * 64
PDF = b"%PDF-1.7\n%%EOF\n"
SUBMISSION_ID = f"{SUBMITTER_ID}_2025-01-01_00000000"


@pytest.fixture
def s3() -> Iterator[Any]:
    """The test config's archives and inbox in a mocked S3, with the submission archived and uploaded."""
    with mock_aws():
        client = boto3.client("s3", region_name="us-east-1")
        for bucket in ("consented", "non_consented", "inbox"):
            client.create_bucket(Bucket=bucket)
        for key in ("files/a.fastq.gz.c4gh", "metadata/metadata.json", "logs/progress_upload.cjson"):
            client.put_object(Bucket="consented", Key=f"{SUBMISSION_ID}/{key}", Body=b"x" * 10)
        for key in ("files/a.fastq.gz.c4gh", "metadata/metadata.json"):
            client.put_object(Bucket="inbox", Key=f"{SUBMISSION_ID}/{key}", Body=b"x" * 10)
        yield client


def _keys(s3: Any, bucket: str) -> set[str]:
    return {o["Key"] for o in s3.list_objects_v2(Bucket=bucket).get("Contents", [])}


def _request(**overrides) -> dict:
    request = {
        "instructed_at": datetime.date(2026, 11, 21),
        "submitter_id": SUBMITTER_ID,
        "requester_name": "Dr. Erika Mustermann",
        "requester_email": "erika.mustermann@example.org",
        "requested_at": datetime.date(2026, 11, 20),
        "reason": "Widerruf der Einwilligung",
        "data_node_id": GRZ_ID,
        "tan_g": [TAN_G_A],
        "comment": "",
    }
    request.update(overrides)
    return request


def _add(
    db: SubmissionDb,
    suffix: str,
    tan_g: str,
    submitter_id: str = SUBMITTER_ID,
    data_node_id: str = GRZ_ID,
    states: tuple[SubmissionStateEnum, ...] = (),
) -> str:
    submission_id = f"{submitter_id}_2025-01-01_{suffix}"
    db.add_submission(submission_id)
    db.modify_submission(submission_id, "tan_g", tan_g)
    db.modify_submission(submission_id, "submitter_id", submitter_id)
    db.modify_submission(submission_id, "data_node_id", data_node_id)
    for state in states:
        db.update_submission_state(submission_id, state)
    return submission_id


def _resolve(db: SubmissionDb, **overrides):
    return resolve_deletion_request(DeletionRequest.model_validate(_request(**overrides)), db, GRZ_ID)


class TestResolve:
    @pytest.mark.parametrize(
        ("states", "outcome"),
        [
            pytest.param((SubmissionStateEnum.REPORTED,), DeletionOutcome.DELETE, id="reported"),
            pytest.param((SubmissionStateEnum.DELETING,), DeletionOutcome.RESUME, id="interrupted"),
            pytest.param(
                (SubmissionStateEnum.DELETING, SubmissionStateEnum.DELETED),
                DeletionOutcome.ALREADY_DELETED,
                id="deleted",
            ),
        ],
    )
    def test_outcome_by_state(self, db: SubmissionDb, states, outcome: DeletionOutcome):
        submission_id = _add(db, "00000000", TAN_G_A, states=states)

        resolved = _resolve(db)

        assert [(r.tan_g, r.submission.id, r.outcome) for r in resolved.tan_gs] == [(TAN_G_A, submission_id, outcome)]
        expected_to_delete = [submission_id] if outcome != DeletionOutcome.ALREADY_DELETED else []
        assert [s.id for s in resolved.to_delete] == expected_to_delete

    def test_unknown_tan_g_is_not_present(self, db: SubmissionDb):
        _add(db, "00000000", TAN_G_A)

        resolved = _resolve(db, tan_g=[TAN_G_A, TAN_G_UNKNOWN])

        assert [r.outcome for r in resolved.tan_gs] == [DeletionOutcome.DELETE, DeletionOutcome.NOT_PRESENT]
        assert resolved.tan_gs[1].submission is None

    def test_tan_g_matches_in_either_case(self, db: SubmissionDb):
        submission_id = _add(db, "00000000", TAN_G_A)

        assert _resolve(db, tan_g=[TAN_G_A.upper()]).tan_gs[0].submission.id == submission_id

    def test_refuses_a_request_for_another_grz(self, db: SubmissionDb):
        _add(db, "00000000", TAN_G_A)

        with pytest.raises(DeletionRequestRefusedError) as e:
            _resolve(db, data_node_id="GRZM00006")
        assert e.value.problems == ["The request is addressed to GRZM00006, but this is GRZK00007."]

    def test_refuses_a_submission_of_another_le(self, db: SubmissionDb):
        _add(db, "00000000", TAN_G_A)
        _add(db, "11111111", TAN_G_B, submitter_id="260914050")

        with pytest.raises(DeletionRequestRefusedError) as e:
            _resolve(db, tan_g=[TAN_G_A, TAN_G_B])
        assert len(e.value.problems) == 1
        assert "was submitted by 260914050" in e.value.problems[0]

    def test_refuses_a_submission_of_another_grz(self, db: SubmissionDb):
        _add(db, "00000000", TAN_G_A, data_node_id="GRZM00006")

        with pytest.raises(DeletionRequestRefusedError) as e:
            _resolve(db)
        assert "names GRZM00006 as its GRZ" in e.value.problems[0]

    @pytest.mark.parametrize(
        ("named", "warned"),
        [pytest.param([TAN_G_A], True, id="followup-not-named"), pytest.param([TAN_G_A, TAN_G_B], False, id="both")],
    )
    def test_warns_about_submissions_of_the_case_the_request_does_not_name(
        self, db: SubmissionDb, named: list[str], warned: bool
    ):
        initial = _add(db, "00000000", TAN_G_A)
        followup = _add(db, "11111111", TAN_G_B)
        db.modify_submission(initial, "submission_type", SubmissionType.initial.value)
        db.modify_submission(followup, "submission_type", SubmissionType.followup.value)
        case_id = db.create_case(SUBMITTER_ID, "patient-1").id
        db.set_submission_case(initial, case_id)
        db.set_submission_case(followup, case_id)

        warnings = _resolve(db, tan_g=named).case_warnings

        if warned:
            assert len(warnings) == 1
            assert initial in warnings[0] and followup in warnings[0]
        else:
            assert warnings == []


def _run(config_path: Path, tmp_path: Path, request: dict, *options: str) -> click.testing.Result:
    data_file = tmp_path / "antrag.yaml"
    data_file.write_text(yaml.safe_dump(request, allow_unicode=True))
    pdf_file = tmp_path / "antrag.pdf"
    pdf_file.write_bytes(PDF)
    return click.testing.CliRunner().invoke(
        grzctl.cli.build_cli(),
        [
            "--config",
            str(config_path),
            "delete",
            "run",
            "--data-file",
            str(data_file),
            "--raw-content",
            str(pdf_file),
            *options,
        ],
    )


class TestRunCommand:
    def test_dry_run_shows_the_plan_without_writing_out_tan_gs(
        self, s3: Any, db: SubmissionDb, migrated_database_config_path: Path, tmp_path: Path
    ):
        submission_id = _add(db, "00000000", TAN_G_A, states=(SubmissionStateEnum.REPORTED,))

        result = _run(migrated_database_config_path, tmp_path, _request(tan_g=[TAN_G_A, TAN_G_UNKNOWN]), "--dry-run")

        assert result.exit_code == 0, result.output
        assert submission_id in result.stdout
        assert "3f9a…0c21" in result.stdout
        assert "not present at GRZK00007" in result.stdout
        assert TAN_G_A not in result.output and TAN_G_UNKNOWN not in result.output
        assert "consented archive" in result.stdout
        assert "2 objects, 20 B" in result.stdout
        assert "Dry run: nothing was deleted." in result.stderr
        assert db.get_submission(submission_id).get_latest_state().state == SubmissionStateEnum.REPORTED
        assert len(_keys(s3, "consented")) == 3

    def test_reports_every_problem_of_the_request(self, migrated_database_config_path: Path, tmp_path: Path):
        result = _run(
            migrated_database_config_path,
            tmp_path,
            _request(tan_g=["<FILL IN tanG>"], submitter_id="12345"),
            "--dry-run",
        )

        assert result.exit_code != 0
        assert "tan_g: still has the '<FILL IN ...>' placeholder" in result.stderr
        assert "submitter_id:" in result.stderr

    def test_refused_request_names_why(self, migrated_database_config_path: Path, tmp_path: Path):
        result = _run(migrated_database_config_path, tmp_path, _request(data_node_id="GRZM00006"), "--dry-run")

        assert result.exit_code != 0
        assert "cannot be carried out" in result.stderr
        assert "addressed to GRZM00006" in result.stderr

    def test_rejects_an_attachment_that_is_no_pdf_or_png(self, migrated_database_config_path: Path, tmp_path: Path):
        data_file = tmp_path / "antrag.yaml"
        data_file.write_text(yaml.safe_dump(_request(), allow_unicode=True))
        text_file = tmp_path / "antrag.txt"
        text_file.write_text("not a pdf")

        result = click.testing.CliRunner().invoke(
            grzctl.cli.build_cli(),
            [
                "--config",
                str(migrated_database_config_path),
                "delete",
                "run",
                "--data-file",
                str(data_file),
                "--raw-content",
                str(text_file),
                "--dry-run",
            ],
        )

        assert result.exit_code != 0


def _states(db: SubmissionDb, submission_id: str) -> list[SubmissionStateEnum]:
    return [state.state for state in sorted(db.get_submission(submission_id).states, key=lambda s: (s.timestamp, s.id))]


class TestRunDeletes:
    def test_deletes_the_submission_and_records_it(
        self, s3: Any, db: SubmissionDb, migrated_database_config_path: Path, tmp_path: Path
    ):
        _add(db, "00000000", TAN_G_A, states=(SubmissionStateEnum.REPORTED,))
        db.modify_submission(SUBMISSION_ID, "local_case_id", "patient-1")

        result = _run(migrated_database_config_path, tmp_path, _request(), "--yes-i-really-mean-it")

        assert result.exit_code == 0, result.output
        assert _keys(s3, "consented") == {f"{SUBMISSION_ID}/logs/progress_upload.cjson", f"{SUBMISSION_ID}/deleted"}
        assert _keys(s3, "non_consented") == set()
        assert _keys(s3, "inbox") == {f"{SUBMISSION_ID}/metadata/metadata.json", f"{SUBMISSION_ID}/cleaned"}

        submission = db.get_submission(SUBMISSION_ID)
        assert _states(db, SUBMISSION_ID)[-2:] == [SubmissionStateEnum.DELETING, SubmissionStateEnum.DELETED]
        assert submission.local_case_id == f"deleted-{SUBMISSION_ID}"
        assert submission.tan_g == TAN_G_A
        deleted = submission.get_latest_state()
        assert [entry["location"] for entry in deleted.data["storage"]] == ["consented archive", "inbox 'inbox'"]
        assert deleted.data["storage"][0]["logs_kept"] == 1

        (change_request,) = db.list_change_requests()[0].changes
        assert change_request.request_raw_content == PDF
        assert change_request.request_raw_content_type == RequestRawContentType.PDF
        assert change_request.data["instructed_at"] == "2026-11-21"
        assert TAN_G_A not in json.dumps(change_request.data) + json.dumps(deleted.data)

    def test_stops_when_object_lock_keeps_objects(
        self, s3: Any, db: SubmissionDb, migrated_database_config_path: Path, tmp_path: Path
    ):
        s3.delete_bucket(Bucket="non_consented")
        s3.create_bucket(Bucket="non_consented", ObjectLockEnabledForBucket=True)
        _add(db, "00000000", TAN_G_A, states=(SubmissionStateEnum.REPORTED,))

        result = _run(migrated_database_config_path, tmp_path, _request(), "--yes-i-really-mean-it")

        assert result.exit_code != 0
        assert "object lock" in result.stderr
        assert len(_keys(s3, "consented")) == 3
        assert _states(db, SUBMISSION_ID) == [SubmissionStateEnum.REPORTED]
        assert db.list_change_requests() == []

    def test_a_failed_deletion_stays_deleting_and_a_second_run_finishes_it(
        self,
        s3: Any,
        db: SubmissionDb,
        migrated_database_config_path: Path,
        tmp_path: Path,
        monkeypatch: pytest.MonkeyPatch,
    ):
        _add(db, "00000000", TAN_G_A, states=(SubmissionStateEnum.REPORTED,))
        with monkeypatch.context() as patch:
            patch.setattr(execute, "delete_from_inbox", lambda *args: None)
            first = _run(migrated_database_config_path, tmp_path, _request(), "--yes-i-really-mean-it")

        assert first.exit_code != 0
        assert "the inbox 'inbox' still holds 2 objects" in first.output
        assert "stays in state Deleting" in first.output
        assert _states(db, SUBMISSION_ID)[-1] == SubmissionStateEnum.DELETING

        second = _run(migrated_database_config_path, tmp_path, _request(), "--yes-i-really-mean-it")

        assert second.exit_code == 0, second.output
        assert _states(db, SUBMISSION_ID)[-1] == SubmissionStateEnum.DELETED
        assert len(db.list_change_requests()[0].changes) == 1

    def test_leaves_submissions_already_deleted_or_not_present(
        self, s3: Any, db: SubmissionDb, migrated_database_config_path: Path, tmp_path: Path
    ):
        _add(db, "00000000", TAN_G_A, states=(SubmissionStateEnum.DELETING, SubmissionStateEnum.DELETED))

        result = _run(
            migrated_database_config_path, tmp_path, _request(tan_g=[TAN_G_A, TAN_G_UNKNOWN]), "--yes-i-really-mean-it"
        )

        assert result.exit_code == 0, result.output
        assert "Nothing to delete." in result.stderr
        assert len(_keys(s3, "consented")) == 3

    def test_answering_no_deletes_nothing(
        self, s3: Any, db: SubmissionDb, migrated_database_config_path: Path, tmp_path: Path
    ):
        _add(db, "00000000", TAN_G_A, states=(SubmissionStateEnum.REPORTED,))
        data_file = tmp_path / "antrag.yaml"
        data_file.write_text(yaml.safe_dump(_request(), allow_unicode=True))
        pdf_file = tmp_path / "antrag.pdf"
        pdf_file.write_bytes(PDF)

        result = click.testing.CliRunner().invoke(
            grzctl.cli.build_cli(),
            [
                "--config",
                str(migrated_database_config_path),
                "delete",
                "run",
                "--data-file",
                str(data_file),
                "--raw-content",
                str(pdf_file),
            ],
            input="n\n",
        )

        assert result.exit_code != 0
        assert len(_keys(s3, "consented")) == 3
        assert _states(db, SUBMISSION_ID) == [SubmissionStateEnum.REPORTED]


def test_a_failed_start_records_no_request(
    s3: Any,
    db: SubmissionDb,
    migrated_database_config_path: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    """The request and the state Deleting are written together, so a rerun cannot record the request twice."""
    _add(db, "00000000", TAN_G_A, states=(SubmissionStateEnum.REPORTED,))
    real_update = SubmissionDb.update_submission_state

    def failing_deleting(self, submission_id, state, *args, **kwargs):
        if state == SubmissionStateEnum.DELETING:
            raise RuntimeError("database gone")
        return real_update(self, submission_id, state, *args, **kwargs)

    with monkeypatch.context() as patch:
        patch.setattr(SubmissionDb, "update_submission_state", failing_deleting)
        failed = _run(migrated_database_config_path, tmp_path, _request(), "--yes-i-really-mean-it")

    assert failed.exit_code != 0
    assert db.list_change_requests() == []

    rerun = _run(migrated_database_config_path, tmp_path, _request(), "--yes-i-really-mean-it")

    assert rerun.exit_code == 0, rerun.output
    assert len(db.list_change_requests()[0].changes) == 1
