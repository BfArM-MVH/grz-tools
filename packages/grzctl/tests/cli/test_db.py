"""
Tests for grzctl db subcommand
"""

import hashlib
import json
import random
from datetime import UTC, date, datetime
from operator import attrgetter
from pathlib import Path
from textwrap import dedent

import click.testing
import grzctl.cli
import pytest
import sqlalchemy
import yaml
from grz_db.models.submission import Submission, SubmissionDb
from grz_pydantic_models.submission.metadata import REDACTED_TAN, GrzSubmissionMetadata
from grzctl.models.config import DbConfig


def test_all_migrations(blank_initial_database_config_path):
    """Database migrations should work all the way from the oldest supported to the latest version."""
    # add some test data
    config = DbConfig.from_path(blank_initial_database_config_path)
    tan_g = "a2b6c3d9e8f7123456789abcdef0123456789abcdef0123456789abcdef01234"
    pseudonym = "CASE12345"
    submission_id = "123456789_2024-11-08_d0f805c5"
    engine = sqlalchemy.create_engine(config.db.database_url)
    with engine.connect() as connection:
        connection.execute(
            sqlalchemy.insert(Submission),
            {"tan_g": tan_g, "pseudonym": pseudonym, "id": submission_id},
        )
        connection.execute(
            sqlalchemy.text(
                """
                INSERT INTO submission_states (state, data, timestamp, submission_id, author_name, signature)
                VALUES (:state, NULL, :timestamp, :submission_id, :author_name, :signature)
                """
            ),
            {
                "state": "QCING",
                "timestamp": datetime.now(UTC),
                "submission_id": submission_id,
                "author_name": "alice",
                "signature": "dummy",
            },
        )
        connection.commit()

    # ensure db command raises appropriate error before migration
    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()
    args_common = ["db", "--config-file", blank_initial_database_config_path]
    result_premature_list = runner.invoke(cli, [*args_common, "list"])
    assert result_premature_list.exit_code != 0
    assert "Database not at latest schema" in result_premature_list.stderr

    # run the migration
    result_upgrade = runner.invoke(cli, [*args_common, "upgrade"])
    assert result_upgrade.exit_code == 0, result_upgrade.stderr

    # check the test data
    result_show = runner.invoke(cli, [*args_common, "submission", "show", submission_id])
    assert result_show.exit_code == 0, result_show.stderr
    # shorter than tanG and less likely to be truncated in various terminal widths
    assert pseudonym in result_show.stdout, result_show.stdout

    db = SubmissionDb(db_url=config.db.database_url, author=None)
    submission = db.get_submission(submission_id)
    assert submission is not None
    assert submission.selected_for_qc is True


def test_populate(blank_database_config_path: Path, test_metadata_path: Path):
    args_common = ["db", "--config-file", blank_database_config_path]
    metadata = GrzSubmissionMetadata.model_validate_json(test_metadata_path.read_text())

    runner = click.testing.CliRunner(catch_exceptions=False)
    cli = grzctl.cli.build_cli()
    result_add = runner.invoke(cli, [*args_common, "submission", "add", metadata.submission_id])
    assert result_add.exit_code == 0, result_add.stderr

    result_populate = runner.invoke(
        cli, [*args_common, "submission", "populate", metadata.submission_id, str(test_metadata_path), "--no-confirm"]
    )
    assert result_populate.exit_code == 0, result_populate.stderr

    result_show = runner.invoke(cli, [*args_common, "submission", "show", metadata.submission_id])
    assert result_show.exit_code == 0, result_show.stderr
    # shorter than tanG and less likely to be truncated in various terminal widths
    assert metadata.submission.local_case_id in result_show.stdout, result_show.stdout

    config = DbConfig.from_path(blank_database_config_path)
    db = SubmissionDb(db_url=config.db.database_url, author=None)

    submission = db.get_submission(metadata.submission_id)
    assert submission.pseudonym == metadata.submission.local_case_id

    # check that the consent records were populated
    meta_father = metadata.donors[1]
    db_father = db.get_donors(submission_id=metadata.submission_id, pseudonym=meta_father.donor_pseudonym)[0]
    assert {
        consent.no_scope_justification for consent in meta_father.research_consents
    } == db_father.research_consent_missing_justifications

    assert submission.submission_metadata is not None
    assert submission.submission_metadata == metadata.to_redacted_dict()

    assert submission.submission_size == metadata.get_submission_size()


def test_populate_date(blank_database_config_path: Path, test_metadata_path: Path):
    db_args = ["db", "--config-file", blank_database_config_path]
    changed_date = date(2026, 1, 1)
    metadata = GrzSubmissionMetadata.model_validate_json(test_metadata_path.read_text())

    runner = click.testing.CliRunner(catch_exceptions=False)
    cli = grzctl.cli.build_cli()
    result_add = runner.invoke(cli, [*db_args, "submission", "add", metadata.submission_id])
    assert result_add.exit_code == 0, result_add.stderr

    result_populate = runner.invoke(
        cli,
        [
            *db_args,
            "submission",
            "populate",
            metadata.submission_id,
            str(test_metadata_path),
            "--no-confirm",
            "--submission_date",
            changed_date.strftime("%Y-%m-%d"),
        ],
    )
    assert result_populate.exit_code == 0, result_populate.stderr

    result_show = runner.invoke(cli, [*db_args, "submission", "show", metadata.submission_id])
    assert result_show.exit_code == 0, result_show.stderr
    # shorter than tanG and less likely to be truncated in various terminal widths
    assert metadata.submission.local_case_id in result_show.stdout, result_show.stdout

    config = DbConfig.from_path(blank_database_config_path)
    db = SubmissionDb(db_url=config.db.database_url, author=None)

    submission = db.get_submission(metadata.submission_id)
    assert submission.submission_date == changed_date


def test_populate_redacted(tmp_path: Path, blank_database_config_path: Path, test_metadata_path: Path):
    args_common = ["db", "--config-file", blank_database_config_path]
    metadata = GrzSubmissionMetadata.model_validate_json(test_metadata_path.read_text())

    # compute submission ID _before_ tanG is redacted (changing the property return value)
    submission_id = metadata.submission_id

    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()
    result_add = runner.invoke(cli, [*args_common, "submission", "add", submission_id])
    assert result_add.exit_code == 0, result_add.stderr

    # redact the tanG
    metadata.submission.tan_g = REDACTED_TAN
    metadata_path = tmp_path / "metadata.json"
    with open(metadata_path, "w") as metadata_file:
        metadata_file.write(metadata.model_dump_json(by_alias=True))

    with pytest.raises(ValueError, match="Refusing to populate a seemingly-redacted TAN"):
        _ = runner.invoke(
            cli,
            [*args_common, "submission", "populate", submission_id, str(metadata_path), "--no-confirm"],
            catch_exceptions=False,
        )


def test_repopulate(blank_database_config_path: Path, tmp_path: Path, test_metadata_path: Path):
    """
    Repopulating a database should work, including when:
    - two donors from different submitters have the same pseudonym.
    - other tricky situations that may be added in the future.
    """
    rng = random.Random()
    rng.seed(42)
    changed_date = date(2026, 1, 1)

    args_common = ["db", "--config-file", blank_database_config_path]
    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()

    metadata_raw = json.loads(test_metadata_path.read_text())

    # first submission
    metadata_raw["submission"]["submitterId"] = "123456789"
    metadata_raw["submission"]["tanG"] = hashlib.sha256(rng.randbytes(128)).hexdigest()
    metadata_s1 = GrzSubmissionMetadata.model_validate_json(json.dumps(metadata_raw))

    result_add_s1 = runner.invoke(cli, [*args_common, "submission", "add", metadata_s1.submission_id])
    assert result_add_s1.exit_code == 0, result_add_s1.stderr

    metadata_s1_dump_path = tmp_path / "metadata.s1.json"
    with open(metadata_s1_dump_path, "w") as metadata_s1_file:
        json.dump(metadata_raw, metadata_s1_file)

    result_populate_s1 = runner.invoke(
        cli,
        [*args_common, "submission", "populate", metadata_s1.submission_id, str(metadata_s1_dump_path), "--no-confirm"],
    )
    assert result_populate_s1.exit_code == 0, result_populate_s1.stderr

    # second submission with same pseudonym from different submitter
    metadata_raw["submission"]["submitterId"] = "987654321"
    metadata_raw["submission"]["tanG"] = hashlib.sha256(rng.randbytes(128)).hexdigest()
    metadata_s2 = GrzSubmissionMetadata.model_validate_json(json.dumps(metadata_raw))

    result_add_s2 = runner.invoke(cli, [*args_common, "submission", "add", metadata_s2.submission_id])
    assert result_add_s2.exit_code == 0, result_add_s2.stderr

    metadata_s2_dump_path = tmp_path / "metadata.s2.json"
    with open(metadata_s2_dump_path, "w") as metadata_s2_file:
        json.dump(metadata_raw, metadata_s2_file)

    result_populate_s2 = runner.invoke(
        cli,
        [*args_common, "submission", "populate", metadata_s2.submission_id, str(metadata_s2_dump_path), "--no-confirm"],
    )
    assert result_populate_s2.exit_code == 0, result_populate_s2.stderr

    # repopulate s1 with redacted metadata and revoked consent
    metadata_raw["submission"]["submitterId"] = "123456789"
    metadata_raw["submission"]["tanG"] = REDACTED_TAN
    metadata_raw["submission"]["localCaseId"] = ""
    metadata_raw["donors"][0]["researchConsents"][0]["scope"] = None
    metadata_raw["donors"][0]["researchConsents"][0]["noScopeJustification"] = "other patient-related reason"

    metadata_s1_mod_dump_path = tmp_path / "metadata.s1.modified.json"
    with open(metadata_s1_mod_dump_path, "w") as metadata_s1_mod_file:
        json.dump(metadata_raw, metadata_s1_mod_file)

    # need to use metadata_s1 submission_id since tanG is now redacted
    result_repopulate_s1 = runner.invoke(
        cli,
        [
            *args_common,
            "submission",
            "populate",
            metadata_s1.submission_id,
            str(metadata_s1_mod_dump_path),
            "--no-confirm",
            "--ignore-field",
            "tan_g",
            "--ignore-field",
            "pseudonym",
            "--submission_date",
            changed_date.strftime("%Y-%m-%d"),
        ],
    )
    assert result_repopulate_s1.exit_code == 0, result_repopulate_s1.stderr

    # sanity check the database
    with open(blank_database_config_path, encoding="utf-8") as blank_database_config_file:
        config = yaml.load(blank_database_config_file, Loader=yaml.Loader)
    db = SubmissionDb(db_url=config["db"]["database_url"], author=None)

    submissions = db.list_submissions(limit=None)
    assert len(submissions) == 2, "Expected two submissions in database"

    donors_s1 = sorted(db.get_donors(metadata_s1.submission_id), key=attrgetter("pseudonym"))
    assert len(donors_s1) == 2, "Expected two donors in submission 1"

    assert donors_s1[1].pseudonym == "index"
    assert not donors_s1[1].research_consented
    assert donors_s1[1].research_consent_missing_justifications == {"other patient-related reason"}

    submission_s1 = db.get_submission(metadata_s1.submission_id)
    assert submission_s1.submission_date == changed_date

    donors_s2 = sorted(db.get_donors(metadata_s2.submission_id), key=attrgetter("pseudonym"))
    assert len(donors_s2) == 2, "Expected two donors in submission 2"


def test_populate_qc(blank_database_config_path: Path, tmp_path: Path, test_metadata_path: Path):
    args_common = ["db", "--config-file", blank_database_config_path]
    metadata = GrzSubmissionMetadata.model_validate_json(test_metadata_path.read_text())

    runner = click.testing.CliRunner(catch_exceptions=False)
    cli = grzctl.cli.build_cli()
    result_add = runner.invoke(cli, [*args_common, "submission", "add", metadata.submission_id])
    assert result_add.exit_code == 0, result_add.stderr

    # populate submission + donors first to satisfy foreign key constraints
    metadata_raw = json.loads(test_metadata_path.read_text())

    metadata_dump_path = tmp_path / "metadata.json"
    with open(metadata_dump_path, "w") as metadata_file:
        json.dump(metadata_raw, metadata_file)

    metadata = GrzSubmissionMetadata.model_validate_json(json.dumps(metadata_raw))
    result_populate = runner.invoke(
        cli,
        [*args_common, "submission", "populate", metadata.submission_id, str(metadata_dump_path), "--no-confirm"],
    )
    assert result_populate.exit_code == 0, result_populate.stderr

    report_csv_path = tmp_path / "report.csv"
    with open(report_csv_path, "w") as report_csv_file:
        report_csv_file.write(
            dedent("""\
            sampleId,donorPseudonym,labDataName,libraryType,sequenceSubtype,genomicStudySubtype,qualityControlStatus,meanDepthOfCoverage,meanDepthOfCoverageProvided,meanDepthOfCoverageRequired,meanDepthOfCoverageDeviation,meanDepthOfCoverageQCStatus,percentBasesAboveQualityThreshold,qualityThreshold,percentBasesAboveQualityThresholdProvided,percentBasesAboveQualityThresholdRequired,percentBasesAboveQualityThresholdDeviation,percentBasesAboveQualityThresholdQCStatus,targetedRegionsAboveMinCoverage,minCoverage,targetedRegionsAboveMinCoverageProvided,targetedRegionsAboveMinCoverageRequired,targetedRegionsAboveMinCoverageDeviation,targetedRegionsAboveMinCoverageQCStatus
            father1_germline0,bbbbbbbb11111111bbbbbbbb11111111bbbbbbbb11111111bbbbbbbb11111111,Blut DNA normal,wes,germline,tumor+germline,FAIL,19.94,30.0,30.0,-33.53333333333333,THRESHOLD NOT MET,90.67913315460233,30,88.0,85,3.0444694938662797,PASS,1.0,20,1.0,0.8,0.0,PASS
            index0_germline0,index,Blut DNA normal,wes,germline,tumor+germline,PASS,49.84,50.0,30.0,-0.3199999999999932,PASS,90.65953529937444,30,88.0,85,3.022199203834591,PASS,1.0,20,1.0,0.8,0.0,PASS
            index0_somatic0,index,Blut DNA Tumor,wes,somatic,tumor+germline,PASS,49.84,50.0,30.0,-0.3199999999999932,PASS,90.65953529937444,30,88.0,85,3.022199203834591,PASS,1.0,20,1.0,0.8,0.0,PASS
            """)
        )

    result_populate = runner.invoke(
        cli,
        [*args_common, "submission", "populate-qc", metadata.submission_id, str(report_csv_path), "--no-confirm"],
    )
    assert result_populate.exit_code == 0, result_populate.stderr

    with open(blank_database_config_path, encoding="utf-8") as blank_database_config_file:
        config = yaml.load(blank_database_config_file, Loader=yaml.Loader)
    db = SubmissionDb(db_url=config["db"]["database_url"], author=None)

    results = db.get_detailed_qc_results(metadata.submission_id)
    assert len(results) == 3
    father_result = next(r for r in results if r.lab_datum_id == "father1_germline0")
    assert not father_result.mean_depth_of_coverage_passed_qc


def test_update_error_confirm(blank_database_config_path: Path, test_metadata_path: Path):
    """Database should confirm before updating a submission from an Error state."""
    args_common = ["db", "--config-file", blank_database_config_path]
    metadata = GrzSubmissionMetadata.model_validate_json(test_metadata_path.read_text())

    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()
    result_add = runner.invoke(cli, [*args_common, "submission", "add", metadata.submission_id])
    assert result_add.exit_code == 0, result_add.stderr

    result_update1 = runner.invoke(cli, [*args_common, "submission", "update", metadata.submission_id, "Error"])
    assert result_update1.exit_code == 0, result_update1.output

    result_update2 = runner.invoke(cli, [*args_common, "submission", "update", metadata.submission_id, "Validated"])
    assert result_update2.exit_code != 0, result_update2.output
    assert (
        "Submission is currently in an 'Error' state. Are you sure you want to set it to 'Validated'?"
        in result_update2.output
    )

    result_update3 = runner.invoke(
        cli, [*args_common, "submission", "update", "--ignore-error-state", metadata.submission_id, "Validated"]
    )
    assert result_update3.exit_code == 0, result_update3.output


def test_list_sort(blank_database_config_path: Path):
    """
    List command should sort in the expected order:
    0. null latest state timestamp and null submission date
    1. latest state timestamp if not null, otherwise submission date
    """
    args_common = ["db", "--config-file", blank_database_config_path]

    expected_ordering = [
        {"id": "123456789_2025-07-01_a1b2c3d4"},
        {"id": "123456789_2025-07-01_a1b2c3d6", "date": "2025-07-01", "add_a_state": True},
        {"id": "123456789_2025-07-01_a1b2c3d5", "date": "2025-07-05"},
    ]

    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()
    for submission in expected_ordering:
        result_add = runner.invoke(cli, [*args_common, "submission", "add", submission["id"]])
        assert result_add.exit_code == 0, result_add.stderr

        if (submission_date := submission.get("date", None)) is not None:
            result_modify = runner.invoke(
                cli, [*args_common, "submission", "modify", submission["id"], "submission_date", submission_date]
            )
            assert result_modify.exit_code == 0, result_modify.stderr

        if submission.get("add_a_state", False):
            result_modify = runner.invoke(cli, [*args_common, "submission", "update", submission["id"], "Uploaded"])
            assert result_modify.exit_code == 0, result_modify.stderr

    result_list = runner.invoke(cli, [*args_common, "list", "--json"])
    assert result_list.exit_code == 0, result_list.stderr

    result_list_parsed = json.loads(result_list.stdout)
    for i, submission in enumerate(expected_ordering):
        assert submission["id"] == result_list_parsed[i]["id"]


def test_submission_show_json(blank_database_config_path: Path, test_metadata_path: Path):
    """
    `grzctl db submission show --json` should return machine-readable JSON
    with submission metadata and state history.
    """
    args_common = ["db", "--config-file", blank_database_config_path]

    metadata = GrzSubmissionMetadata.model_validate_json(test_metadata_path.read_text())

    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()

    # add submission
    result_add = runner.invoke(cli, [*args_common, "submission", "add", metadata.submission_id])
    assert result_add.exit_code == 0, result_add.stderr

    # populate submission
    result_populate = runner.invoke(
        cli,
        [*args_common, "submission", "populate", metadata.submission_id, str(test_metadata_path), "--no-confirm"],
    )
    assert result_populate.exit_code == 0, result_populate.stderr

    # show submission as JSON
    result_show_json = runner.invoke(cli, [*args_common, "submission", "show", "--json", metadata.submission_id])
    assert result_show_json.exit_code == 0, result_show_json.stderr

    parsed = json.loads(result_show_json.stdout)

    # structure checks
    assert parsed == {
        "id": metadata.submission_id,
        "tan_g": metadata.submission.tan_g,
        "pseudonym": metadata.submission.local_case_id,
        "submission_date": metadata.submission.submission_date.isoformat()
        if metadata.submission.submission_date
        else None,
        "submission_size": metadata.get_submission_size(),
        "submission_type": metadata.submission.submission_type,
        "submission_metadata": metadata.to_redacted_dict(),
        "submitter_id": metadata.submission.submitter_id,
        "data_node_id": metadata.submission.genomic_data_center_id,
        "coverage_type": metadata.submission.coverage_type,
        "disease_type": metadata.submission.disease_type,
        "basic_qc_passed": None,
        "consented": metadata.consents_to_research(date=date.today()),
        "selected_for_qc": None,
        "detailed_qc_passed": None,
        "genomic_study_type": metadata.submission.genomic_study_type,
        "genomic_study_subtype": metadata.submission.genomic_study_subtype,
        "states": [],
    }


def test_list_filter_modes_and_multiple_states(blank_database_config_path: Path):
    args_common = ["db", "--config-file", blank_database_config_path]

    sub_latest_error = "123456789_2025-07-01_b1b2c3d1"
    sub_latest_qced = "123456789_2025-07-01_b1b2c3d2"
    sub_latest_downloaded = "123456789_2025-07-01_b1b2c3d3"
    sub_history_error_latest_uploaded = "123456789_2025-07-01_b1b2c3d4"
    sub_no_state = "123456789_2025-07-01_b1b2c3d5"

    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()

    for submission_id in [
        sub_latest_error,
        sub_latest_qced,
        sub_latest_downloaded,
        sub_history_error_latest_uploaded,
        sub_no_state,
    ]:
        result_add = runner.invoke(cli, [*args_common, "submission", "add", submission_id])
        assert result_add.exit_code == 0, result_add.stderr

    update_invocations = [
        [*args_common, "submission", "update", sub_latest_error, "Error"],
        [*args_common, "submission", "update", sub_latest_qced, "QCed"],
        [*args_common, "submission", "update", sub_latest_downloaded, "Downloaded"],
        [*args_common, "submission", "update", sub_history_error_latest_uploaded, "Uploaded"],
        [*args_common, "submission", "update", sub_history_error_latest_uploaded, "Error"],
    ]
    for invoke_args in update_invocations:
        result_update = runner.invoke(cli, invoke_args)
        assert result_update.exit_code == 0, result_update.stderr

    # Special transition from Error -> Uploaded, explicitly allowed by flag.
    result_ignore_error = runner.invoke(
        cli,
        [
            *args_common,
            "submission",
            "update",
            "--ignore-error-state",
            sub_history_error_latest_uploaded,
            "Uploaded",
        ],
    )
    assert result_ignore_error.exit_code == 0, result_ignore_error.stderr

    result_all = runner.invoke(cli, [*args_common, "list", "--json"])
    assert result_all.exit_code == 0, result_all.stderr
    parsed_all = json.loads(result_all.stdout)
    assert {item["id"] for item in parsed_all} == {
        sub_latest_error,
        sub_latest_qced,
        sub_latest_downloaded,
        sub_history_error_latest_uploaded,
        sub_no_state,
    }

    # default mode is latest
    result_error_latest_default = runner.invoke(cli, [*args_common, "list", "--json", "--state", "error"])
    assert result_error_latest_default.exit_code == 0, result_error_latest_default.stderr
    parsed_error_latest_default = json.loads(result_error_latest_default.stdout)
    assert {item["id"] for item in parsed_error_latest_default} == {sub_latest_error}
    assert parsed_error_latest_default[0]["latest_state"]["state"] == "Error"

    # any mode includes submissions that had Error in history but not as latest state
    result_error_any = runner.invoke(
        cli,
        [*args_common, "list", "--json", "--state", "error", "--filter-mode", "any"],
    )
    assert result_error_any.exit_code == 0, result_error_any.stderr
    parsed_error_any = json.loads(result_error_any.stdout)
    assert {item["id"] for item in parsed_error_any} == {sub_latest_error, sub_history_error_latest_uploaded}
    history_match = next(filter(lambda item: item["id"] == sub_history_error_latest_uploaded, parsed_error_any))
    assert history_match["latest_state"]["state"] == "Uploaded"

    # filter mode parsing remains case-insensitive
    result_error_any_upper = runner.invoke(
        cli,
        [*args_common, "list", "--json", "--state", "error", "--filter-mode", "ANY"],
    )
    assert result_error_any_upper.exit_code == 0, result_error_any_upper.stderr
    parsed_error_any_upper = json.loads(result_error_any_upper.stdout)
    assert {item["id"] for item in parsed_error_any_upper} == {sub_latest_error, sub_history_error_latest_uploaded}

    # multiple filters in default latest mode are OR-ed on latest state
    result_multi_latest = runner.invoke(
        cli,
        [*args_common, "list", "--json", "--state", "error", "--state", "qced"],
    )
    assert result_multi_latest.exit_code == 0, result_multi_latest.stderr
    parsed_multi_latest = json.loads(result_multi_latest.stdout)
    assert {item["id"] for item in parsed_multi_latest} == {sub_latest_error, sub_latest_qced}

    # multiple filters in any mode are OR-ed across all historic states
    result_multi_any = runner.invoke(
        cli,
        [*args_common, "list", "--json", "--state", "error", "--state", "qced", "--filter-mode", "any"],
    )
    assert result_multi_any.exit_code == 0, result_multi_any.stderr
    parsed_multi_any = json.loads(result_multi_any.stdout)
    assert {item["id"] for item in parsed_multi_any} == {
        sub_latest_error,
        sub_history_error_latest_uploaded,
        sub_latest_qced,
    }

    # "--state" supports repeated values
    result_state_alias = runner.invoke(
        cli,
        [*args_common, "list", "--json", "--state", "downloaded"],
    )
    assert result_state_alias.exit_code == 0, result_state_alias.stderr
    parsed_state_alias = json.loads(result_state_alias.stdout)
    assert {item["id"] for item in parsed_state_alias} == {sub_latest_downloaded}


_DELETE_CHANGE_REQUEST_DATA = {
    "requester_name": "Erika Mustermann",
    "requester_email": "demo.requester@example.org",
    "requested_at": "2026-03-27",
    "request_email_content": "Liebe Kolleginnen,\nbitte löschen.\nVielen Dank,\nErika",
}


def _assert_audit_columns_match(change_log, expected: dict) -> None:
    """Check the audit columns on a ChangeRequestLog row match the expected input dict."""
    assert change_log.requester_name == expected["requester_name"]
    assert change_log.requester_email == expected["requester_email"]
    assert change_log.requested_at.isoformat() == expected["requested_at"]
    assert change_log.request_email_content == expected["request_email_content"]


def _add_submission(runner: click.testing.CliRunner, cli, args_common: list, submission_id: str) -> None:
    result = runner.invoke(cli, [*args_common, "submission", "add", submission_id])
    assert result.exit_code == 0, result.stderr


def test_change_request_delete_succeeds_with_yaml_data_file(blank_database_config_path: Path, tmp_path: Path):
    args_common = ["db", "--config-file", blank_database_config_path]
    submission_id = "260840108_2025-12-16_cc9973f0"
    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()
    _add_submission(runner, cli, args_common, submission_id)

    data_file = tmp_path / "delete.yaml"
    data_file.write_text(yaml.safe_dump(_DELETE_CHANGE_REQUEST_DATA, allow_unicode=True))

    result = runner.invoke(
        cli,
        [
            *args_common,
            "submission",
            "change-request",
            submission_id,
            "Delete",
            "--data-file",
            str(data_file),
        ],
    )
    assert result.exit_code == 0, result.stderr
    assert "has undergone a change request" in result.stderr.replace("\n", " ")
    assert "Delete" in result.stderr

    config = DbConfig.from_path(blank_database_config_path)
    db = SubmissionDb(db_url=config.db.database_url, author=None)
    submissions_with_changes = [s for s in db.list_change_requests() if s.id == submission_id]
    assert len(submissions_with_changes) == 1
    changes = submissions_with_changes[0].changes
    assert len(changes) == 1
    _assert_audit_columns_match(changes[0], _DELETE_CHANGE_REQUEST_DATA)
    assert changes[0].data is None


def test_change_request_delete_succeeds_with_inline_data(blank_database_config_path: Path):
    """`--data` accepts valid Delete data and stores it the same way `--data-file` does."""
    args_common = ["db", "--config-file", blank_database_config_path]
    submission_id = "260840108_2025-12-16_cc9973f0"
    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()
    _add_submission(runner, cli, args_common, submission_id)

    result = runner.invoke(
        cli,
        [
            *args_common,
            "submission",
            "change-request",
            submission_id,
            "Delete",
            "--data",
            json.dumps(_DELETE_CHANGE_REQUEST_DATA),
        ],
    )
    assert result.exit_code == 0, result.stderr

    config = DbConfig.from_path(blank_database_config_path)
    db = SubmissionDb(db_url=config.db.database_url, author=None)
    submissions_with_changes = [s for s in db.list_change_requests() if s.id == submission_id]
    _assert_audit_columns_match(submissions_with_changes[0].changes[0], _DELETE_CHANGE_REQUEST_DATA)


def test_change_request_delete_inline_data_rejects_placeholder(blank_database_config_path: Path):
    """`--data` is subject to the same placeholder rejection as `--data-file`."""
    args_common = ["db", "--config-file", blank_database_config_path]
    submission_id = "260840108_2025-12-16_cc9973f0"
    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()
    _add_submission(runner, cli, args_common, submission_id)

    bad = {**_DELETE_CHANGE_REQUEST_DATA, "requester_name": "<FILL IN requester name>"}
    result = runner.invoke(
        cli,
        [
            *args_common,
            "submission",
            "change-request",
            submission_id,
            "Delete",
            "--data",
            json.dumps(bad),
        ],
    )
    assert result.exit_code != 0
    assert "FILL IN" in result.stderr or "placeholder" in result.stderr


def test_change_request_delete_inline_data_rejects_missing_field(blank_database_config_path: Path):
    """`--data` is subject to the same required-field check as `--data-file`."""
    args_common = ["db", "--config-file", blank_database_config_path]
    submission_id = "260840108_2025-12-16_cc9973f0"
    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()
    _add_submission(runner, cli, args_common, submission_id)

    incomplete = {k: v for k, v in _DELETE_CHANGE_REQUEST_DATA.items() if k != "requester_email"}
    result = runner.invoke(
        cli,
        [
            *args_common,
            "submission",
            "change-request",
            submission_id,
            "Delete",
            "--data",
            json.dumps(incomplete),
        ],
    )
    assert result.exit_code != 0
    assert "requester_email" in result.stderr


def test_change_request_delete_fails_without_data(blank_database_config_path: Path):
    args_common = ["db", "--config-file", blank_database_config_path]
    submission_id = "260840108_2025-12-16_cc9973f0"
    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()
    _add_submission(runner, cli, args_common, submission_id)

    result = runner.invoke(cli, [*args_common, "submission", "change-request", submission_id, "Delete"])
    assert result.exit_code != 0
    assert "requires audit data" in result.stderr
    assert "requester_name" in result.stderr
    assert "request_email_content" in result.stderr


def test_change_request_delete_fails_with_missing_field(blank_database_config_path: Path, tmp_path: Path):
    args_common = ["db", "--config-file", blank_database_config_path]
    submission_id = "260840108_2025-12-16_cc9973f0"
    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()
    _add_submission(runner, cli, args_common, submission_id)

    incomplete = {k: v for k, v in _DELETE_CHANGE_REQUEST_DATA.items() if k != "requester_email"}
    data_file = tmp_path / "delete.yaml"
    data_file.write_text(yaml.safe_dump(incomplete, allow_unicode=True))

    result = runner.invoke(
        cli,
        [
            *args_common,
            "submission",
            "change-request",
            submission_id,
            "Delete",
            "--data-file",
            str(data_file),
        ],
    )
    assert result.exit_code != 0
    assert "requester_email" in result.stderr


def test_change_request_delete_accepts_and_persists_extra_field(blank_database_config_path: Path, tmp_path: Path):
    """Extras are allowed alongside the required schema fields and stored verbatim."""
    args_common = ["db", "--config-file", blank_database_config_path]
    submission_id = "260840108_2025-12-16_cc9973f0"
    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()
    _add_submission(runner, cli, args_common, submission_id)

    data_with_extra = {**_DELETE_CHANGE_REQUEST_DATA, "internal_note": "received fax 2026-04-30"}
    data_file = tmp_path / "delete.yaml"
    data_file.write_text(yaml.safe_dump(data_with_extra, allow_unicode=True))

    result = runner.invoke(
        cli,
        [
            *args_common,
            "submission",
            "change-request",
            submission_id,
            "Delete",
            "--data-file",
            str(data_file),
        ],
    )
    assert result.exit_code == 0, result.stderr

    config = DbConfig.from_path(blank_database_config_path)
    db = SubmissionDb(db_url=config.db.database_url, author=None)
    submissions_with_changes = [s for s in db.list_change_requests() if s.id == submission_id]
    persisted = submissions_with_changes[0].changes[0]
    _assert_audit_columns_match(persisted, _DELETE_CHANGE_REQUEST_DATA)
    assert persisted.data == {"internal_note": "received fax 2026-04-30"}


def test_change_request_data_and_data_file_mutually_exclusive(blank_database_config_path: Path, tmp_path: Path):
    args_common = ["db", "--config-file", blank_database_config_path]
    submission_id = "260840108_2025-12-16_cc9973f0"
    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()
    _add_submission(runner, cli, args_common, submission_id)

    data_file = tmp_path / "delete.yaml"
    data_file.write_text(yaml.safe_dump(_DELETE_CHANGE_REQUEST_DATA, allow_unicode=True))

    result = runner.invoke(
        cli,
        [
            *args_common,
            "submission",
            "change-request",
            submission_id,
            "Delete",
            "--data",
            json.dumps(_DELETE_CHANGE_REQUEST_DATA),
            "--data-file",
            str(data_file),
        ],
    )
    assert result.exit_code != 0
    assert "mutually exclusive" in result.stderr


def test_change_request_template_subcommand_prints_template():
    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()
    result = runner.invoke(cli, ["change-request-template", "Delete"])
    assert result.exit_code == 0, result.stderr
    expected_audit_fields = {
        "requester_name",
        "requester_email",
        "requested_at",
        "request_email_content",
    }
    for field_name in expected_audit_fields:
        assert field_name in result.stdout
    # template should be valid YAML and include the optional `data` extras section
    parsed = yaml.safe_load(result.stdout)
    assert isinstance(parsed, dict)
    assert expected_audit_fields <= set(parsed.keys())
    assert "data" in parsed


def test_unmodified_template_fails_validation(blank_database_config_path: Path, tmp_path: Path):
    """A user who saves the template and submits it unchanged must not pass validation."""
    args_common = ["db", "--config-file", blank_database_config_path]
    submission_id = "260840108_2025-12-16_cc9973f0"
    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()
    _add_submission(runner, cli, args_common, submission_id)

    template_result = runner.invoke(cli, ["change-request-template", "Delete"])
    assert template_result.exit_code == 0, template_result.stderr

    data_file = tmp_path / "delete.yaml"
    data_file.write_text(template_result.stdout)

    result = runner.invoke(
        cli,
        [
            *args_common,
            "submission",
            "change-request",
            submission_id,
            "Delete",
            "--data-file",
            str(data_file),
        ],
    )
    assert result.exit_code != 0
    # `requested_at: YYYY-MM-DD` fails date coercion before the model validator runs.
    assert "requested_at" in result.stderr


def test_template_with_only_date_filled_in_still_fails(blank_database_config_path: Path, tmp_path: Path):
    """If the user fills in only the date but leaves text placeholders, validation must still fail."""
    args_common = ["db", "--config-file", blank_database_config_path]
    submission_id = "260840108_2025-12-16_cc9973f0"
    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()
    _add_submission(runner, cli, args_common, submission_id)

    template_result = runner.invoke(cli, ["change-request-template", "Delete"])
    partially_filled = template_result.stdout.replace("YYYY-MM-DD", "2026-03-27")
    data_file = tmp_path / "delete.yaml"
    data_file.write_text(partially_filled)

    result = runner.invoke(
        cli,
        [
            *args_common,
            "submission",
            "change-request",
            submission_id,
            "Delete",
            "--data-file",
            str(data_file),
        ],
    )
    assert result.exit_code != 0
    assert "FILL IN" in result.stderr or "placeholder" in result.stderr


def test_change_request_template_for_other_change_types_includes_audit_fields():
    """Audit fields are universal — every change type prints the same scaffold (with type-specific guidance)."""
    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()
    result = runner.invoke(cli, ["change-request-template", "Modify"])
    assert result.exit_code == 0, result.stderr
    for field_name in ("requester_name", "requester_email", "requested_at", "request_email_content"):
        assert field_name in result.stdout
    assert "describe the field/value to modify" in result.stdout


def test_change_request_dry_run_does_not_write(blank_database_config_path: Path, tmp_path: Path):
    args_common = ["db", "--config-file", blank_database_config_path]
    submission_id = "260840108_2025-12-16_cc9973f0"
    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()
    _add_submission(runner, cli, args_common, submission_id)

    data_file = tmp_path / "delete.yaml"
    data_file.write_text(yaml.safe_dump(_DELETE_CHANGE_REQUEST_DATA, allow_unicode=True))

    result = runner.invoke(
        cli,
        [
            *args_common,
            "submission",
            "change-request",
            submission_id,
            "Delete",
            "--data-file",
            str(data_file),
            "--dry-run",
        ],
    )
    assert result.exit_code == 0, result.stderr
    assert "Dry run" in result.stderr
    assert "would register change request" in result.stderr.replace("\n", " ")
    assert "Validated fields" in result.stderr

    config = DbConfig.from_path(blank_database_config_path)
    db = SubmissionDb(db_url=config.db.database_url, author=None)
    submissions_with_changes = list(db.list_change_requests())
    assert submissions_with_changes == []


def test_change_request_dry_run_aborts_when_submission_missing(blank_database_config_path: Path, tmp_path: Path):
    args_common = ["db", "--config-file", blank_database_config_path]
    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()
    data_file = tmp_path / "delete.yaml"
    data_file.write_text(yaml.safe_dump(_DELETE_CHANGE_REQUEST_DATA, allow_unicode=True))

    result = runner.invoke(
        cli,
        [
            *args_common,
            "submission",
            "change-request",
            "DOES_NOT_EXIST",
            "Delete",
            "--data-file",
            str(data_file),
            "--dry-run",
        ],
    )
    assert result.exit_code != 0
    assert "not found" in result.stderr


def test_change_request_dry_run_validates_before_db_check(blank_database_config_path: Path, tmp_path: Path):
    """Schema validation runs first; missing fields fail even with --dry-run."""
    args_common = ["db", "--config-file", blank_database_config_path]
    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()
    incomplete = {k: v for k, v in _DELETE_CHANGE_REQUEST_DATA.items() if k != "requester_email"}
    data_file = tmp_path / "delete.yaml"
    data_file.write_text(yaml.safe_dump(incomplete, allow_unicode=True))

    result = runner.invoke(
        cli,
        [
            *args_common,
            "submission",
            "change-request",
            "DOES_NOT_EXIST",
            "Delete",
            "--data-file",
            str(data_file),
            "--dry-run",
        ],
    )
    assert result.exit_code != 0
    assert "requester_email" in result.stderr


def test_change_request_modify_requires_audit_fields_too(blank_database_config_path: Path, tmp_path: Path):
    """Audit fields are universal — Modify also requires them via --data/--data-file."""
    args_common = ["db", "--config-file", blank_database_config_path]
    submission_id = "260840108_2025-12-16_cc9973f0"
    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()
    _add_submission(runner, cli, args_common, submission_id)

    bare = runner.invoke(cli, [*args_common, "submission", "change-request", submission_id, "Modify"])
    assert bare.exit_code != 0
    assert "requires audit data" in bare.stderr

    data_file = tmp_path / "modify.yaml"
    data_file.write_text(
        yaml.safe_dump(
            {**_DELETE_CHANGE_REQUEST_DATA, "request_email_content": "please modify field X to Y"},
            allow_unicode=True,
        )
    )
    ok = runner.invoke(
        cli,
        [*args_common, "submission", "change-request", submission_id, "Modify", "--data-file", str(data_file)],
    )
    assert ok.exit_code == 0, ok.stderr
    assert "Modify" in ok.stderr


def test_change_request_with_raw_content_pdf(blank_database_config_path: Path, tmp_path: Path):
    """--raw-content provides the binary blob; type is inferred from .pdf extension."""
    args_common = ["db", "--config-file", blank_database_config_path]
    submission_id = "260840108_2025-12-16_cc9973f0"
    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()
    _add_submission(runner, cli, args_common, submission_id)

    pdf_bytes = b"%PDF-1.4 fake pdf body \xde\xad\xbe\xef"
    pdf_path = tmp_path / "request.pdf"
    pdf_path.write_bytes(pdf_bytes)

    audit_only = {k: v for k, v in _DELETE_CHANGE_REQUEST_DATA.items() if k != "request_email_content"}
    data_file = tmp_path / "delete.yaml"
    data_file.write_text(yaml.safe_dump(audit_only, allow_unicode=True))

    result = runner.invoke(
        cli,
        [
            *args_common,
            "submission",
            "change-request",
            submission_id,
            "Delete",
            "--data-file",
            str(data_file),
            "--raw-content",
            str(pdf_path),
        ],
    )
    assert result.exit_code == 0, result.stderr

    config = DbConfig.from_path(blank_database_config_path)
    db = SubmissionDb(db_url=config.db.database_url, author=None)
    persisted = next(s for s in db.list_change_requests() if s.id == submission_id).changes[0]
    assert persisted.request_email_content is None
    assert persisted.request_raw_content == pdf_bytes
    assert persisted.request_raw_content_type.value == "PDF"


def test_change_request_raw_content_unknown_extension_fails(blank_database_config_path: Path, tmp_path: Path):
    """An unknown extension must fail with a clear error since type is only inferred."""
    args_common = ["db", "--config-file", blank_database_config_path]
    submission_id = "260840108_2025-12-16_cc9973f0"
    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()
    _add_submission(runner, cli, args_common, submission_id)

    blob_path = tmp_path / "request.bin"
    blob_path.write_bytes(b"\x00\x01\x02")
    data_file = tmp_path / "delete.yaml"
    data_file.write_text(yaml.safe_dump(_DELETE_CHANGE_REQUEST_DATA, allow_unicode=True))

    result = runner.invoke(
        cli,
        [
            *args_common,
            "submission",
            "change-request",
            submission_id,
            "Delete",
            "--data-file",
            str(data_file),
            "--raw-content",
            str(blob_path),
        ],
    )
    assert result.exit_code != 0
    assert "infer raw-content type" in result.stderr


def test_change_request_raw_content_magic_byte_mismatch_fails(blank_database_config_path: Path, tmp_path: Path):
    """A file with a .pdf extension but non-PDF bytes must fail the magic-byte validator."""
    args_common = ["db", "--config-file", blank_database_config_path]
    submission_id = "260840108_2025-12-16_cc9973f0"
    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()
    _add_submission(runner, cli, args_common, submission_id)

    fake_pdf = tmp_path / "request.pdf"
    fake_pdf.write_bytes(b"not a pdf at all")
    data_file = tmp_path / "delete.yaml"
    data_file.write_text(yaml.safe_dump(_DELETE_CHANGE_REQUEST_DATA, allow_unicode=True))

    result = runner.invoke(
        cli,
        [
            *args_common,
            "submission",
            "change-request",
            submission_id,
            "Delete",
            "--data-file",
            str(data_file),
            "--raw-content",
            str(fake_pdf),
        ],
    )
    assert result.exit_code != 0
    assert "magic bytes" in result.stderr
