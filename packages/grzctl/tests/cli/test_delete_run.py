import datetime
from pathlib import Path

import click.testing
import grzctl.cli
import pytest
import yaml
from grz_db.models.submission import SubmissionDb, SubmissionStateEnum, SubmissionType
from grzctl.commands.delete.resolve import (
    DeletionOutcome,
    DeletionRequestRefusedError,
    resolve_deletion_request,
)
from grzctl.models.deletion_request import DeletionRequest

GRZ_ID = "GRZK00007"
SUBMITTER_ID = "000000000"
TAN_G_A = "3f9a" + "0" * 56 + "0c21"
TAN_G_B = "8b1e" + "1" * 56 + "00d4"
TAN_G_UNKNOWN = "c" * 64
PDF = b"%PDF-1.7\n%%EOF\n"


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
        self, db: SubmissionDb, migrated_database_config_path: Path, tmp_path: Path
    ):
        submission_id = _add(db, "00000000", TAN_G_A, states=(SubmissionStateEnum.REPORTED,))

        result = _run(migrated_database_config_path, tmp_path, _request(tan_g=[TAN_G_A, TAN_G_UNKNOWN]), "--dry-run")

        assert result.exit_code == 0, result.output
        assert submission_id in result.stdout
        assert "3f9a…0c21" in result.stdout
        assert "not present at GRZK00007" in result.stdout
        assert TAN_G_A not in result.output and TAN_G_UNKNOWN not in result.output
        assert "Dry run: nothing was deleted." in result.stderr
        assert db.get_submission(submission_id).get_latest_state().state == SubmissionStateEnum.REPORTED

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
