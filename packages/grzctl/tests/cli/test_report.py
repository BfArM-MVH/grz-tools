import csv
import importlib.resources
import json
import sqlite3
from operator import itemgetter
from pathlib import Path
from textwrap import dedent

import grzctl.cli
from click.testing import CliRunner
from grz_pydantic_models.submission.metadata import GrzSubmissionMetadata
from grzctl.models.config import DbConfig

from .. import resources as test_resources


def test_quarterly_empty(blank_database_config_path: Path, tmp_path: Path):
    """Quarterly reports should work on an empty database."""
    env = {"GRZ_DB__AUTHOR__PRIVATE_KEY_PASSPHRASE": "test"}

    runner = CliRunner(env=env)
    cli = grzctl.cli.build_cli()

    with runner.isolated_filesystem(temp_dir=tmp_path) as report_tmp_dir:
        result_report = runner.invoke(
            cli, ["report", "--config-file", blank_database_config_path, "quarterly"], catch_exceptions=False
        )
        assert result_report.exit_code == 0, result_report.output

    assert (Path(report_tmp_dir) / "Gesamtübersicht.tsv").exists()
    assert (Path(report_tmp_dir) / "Infos_zu_Datensätzen.tsv").exists()
    assert (Path(report_tmp_dir) / "Detailprüfung.tsv").exists()


def test_quarterly(blank_database_config_path: Path, tmp_path: Path):
    """Small test case with a few submissions for quarterly reports."""
    env = {"GRZ_DB__AUTHOR__PRIVATE_KEY_PASSPHRASE": "test"}

    runner = CliRunner(env=env)
    cli = grzctl.cli.build_cli()

    # add and populate first submission to database
    s1_metadata_raw = json.loads((importlib.resources.files(test_resources) / "metadata.json").read_text())
    s1_metadata = GrzSubmissionMetadata.model_validate(s1_metadata_raw)
    s1_metadata_raw["submission"]["submissionType"] = "initial"
    result_add1 = runner.invoke(
        cli, ["db", "--config-file", blank_database_config_path, "submission", "add", s1_metadata.submission_id]
    )
    assert result_add1.exit_code == 0, result_add1.output
    s1_metadata_path = tmp_path / "submission1.metadata.json"
    with open(s1_metadata_path, mode="w", encoding="utf-8") as s1_metadata_file:
        json.dump(s1_metadata_raw, s1_metadata_file)
    result_populate1 = runner.invoke(
        cli,
        [
            "db",
            "--config-file",
            blank_database_config_path,
            "submission",
            "populate",
            "--no-confirm",
            s1_metadata.submission_id,
            str(s1_metadata_path),
        ],
    )
    assert result_populate1.exit_code == 0, result_populate1.output
    result_modify1 = runner.invoke(
        cli,
        [
            "db",
            "--config-file",
            blank_database_config_path,
            "submission",
            "modify",
            s1_metadata.submission_id,
            "basic_qc_passed",
            "yes",
        ],
    )
    assert result_modify1.exit_code == 0, result_modify1.output

    # add a single test submission from another submitter that fails detailed QC
    s2_metadata_raw = json.loads((importlib.resources.files(test_resources) / "metadata.json").read_text())
    s2_metadata_raw["submission"]["submitterId"] = "987654321"
    s2_metadata_raw["submission"]["genomicStudyType"] = "single"
    s2_metadata_raw["submission"]["tanG"] = "d92f44b998916af883c7d4df8a94f06a9f3bf66e3e1c94753a5c310043b2cb24"
    del s2_metadata_raw["donors"][-1]
    s2_metadata_raw["donors"][0]["researchConsents"][0]["scope"] = None
    s2_metadata_raw["donors"][0]["researchConsents"][0]["noScopeJustification"] = "patient refuses to sign consent"
    s2_metadata = GrzSubmissionMetadata.model_validate(s2_metadata_raw)
    result_add2 = runner.invoke(
        cli, ["db", "--config-file", blank_database_config_path, "submission", "add", s2_metadata.submission_id]
    )
    assert result_add2.exit_code == 0, result_add2.output
    s2_metadata_path = tmp_path / "submission2.metadata.json"
    with open(s2_metadata_path, mode="w", encoding="utf-8") as s2_metadata_file:
        json.dump(s2_metadata_raw, s2_metadata_file)
    result_populate2 = runner.invoke(
        cli,
        [
            "db",
            "--config-file",
            blank_database_config_path,
            "submission",
            "populate",
            "--no-confirm",
            s2_metadata.submission_id,
            str(s2_metadata_path),
        ],
    )
    assert result_populate2.exit_code == 0, result_populate2.output
    result_modify2 = runner.invoke(
        cli,
        [
            "db",
            "--config-file",
            blank_database_config_path,
            "submission",
            "modify",
            s2_metadata.submission_id,
            "basic_qc_passed",
            "yes",
        ],
    )
    assert result_modify2.exit_code == 0, result_modify2.output
    result_modify2_2 = runner.invoke(
        cli,
        [
            "db",
            "--config-file",
            blank_database_config_path,
            "submission",
            "modify",
            s2_metadata.submission_id,
            "detailed_qc_passed",
            "no",
        ],
    )
    assert result_modify2_2.exit_code == 0, result_modify2_2.output
    report_csv_path = tmp_path / "submission2.report.csv"
    with open(report_csv_path, "w") as report_csv_file:
        report_csv_file.write(
            dedent("""\
            sampleId,donorPseudonym,labDataName,libraryType,sequenceSubtype,genomicStudySubtype,qualityControlStatus,meanDepthOfCoverage,meanDepthOfCoverageProvided,meanDepthOfCoverageRequired,meanDepthOfCoverageDeviation,meanDepthOfCoverageQCStatus,percentBasesAboveQualityThreshold,qualityThreshold,percentBasesAboveQualityThresholdProvided,percentBasesAboveQualityThresholdRequired,percentBasesAboveQualityThresholdDeviation,percentBasesAboveQualityThresholdQCStatus,targetedRegionsAboveMinCoverage,minCoverage,targetedRegionsAboveMinCoverageProvided,targetedRegionsAboveMinCoverageRequired,targetedRegionsAboveMinCoverageDeviation,targetedRegionsAboveMinCoverageQCStatus
            index0_germline0,index,Blood DNA normal,wes,germline,tumor+germline,PASS,45,50.0,30.0,-10,TOO LOW,90.65953529937444,30,88.0,85,3.022199203834591,PASS,1.0,20,1.0,0.8,0.0,PASS
            index0_somatic0,index,Blood DNA tumor,wes,somatic,tumor+germline,FAIL,49.84,50.0,30.0,-0.3199999999999932,PASS,90.65953529937444,30,88.0,85,3.022199203834591,PASS,1.0,20,1.0,0.8,0.0,PASS
            """)
        )
    result_qc_populate2 = runner.invoke(
        cli,
        [
            "db",
            "--config-file",
            blank_database_config_path,
            "submission",
            "populate-qc",
            s2_metadata.submission_id,
            str(report_csv_path),
            "--no-confirm",
        ],
    )
    assert result_qc_populate2.exit_code == 0, result_qc_populate2.output

    # add correction submission that revokes consent of index patient in first submission, and add deletion change request for original submission
    s3_metadata_raw = json.loads((importlib.resources.files(test_resources) / "metadata.json").read_text())
    s3_metadata_raw["submission"]["submissionType"] = "correction"
    s3_metadata_raw["submission"]["tanG"] = "e8bd8d543a8590d9baf7302dad693ecd77fe12a8760f92ce7be4dddb15681788"
    s3_metadata_raw["donors"][0]["researchConsents"][0]["scope"]["provision"]["provision"] = []
    s3_metadata = GrzSubmissionMetadata.model_validate(s3_metadata_raw)
    result_add3 = runner.invoke(
        cli, ["db", "--config-file", blank_database_config_path, "submission", "add", s3_metadata.submission_id]
    )
    assert result_add3.exit_code == 0, result_add3.output
    s3_metadata_path = tmp_path / "submission3.metadata.json"
    with open(s3_metadata_path, mode="w", encoding="utf-8") as s3_metadata_file:
        json.dump(s3_metadata_raw, s3_metadata_file)
    result_populate3 = runner.invoke(
        cli,
        [
            "db",
            "--config-file",
            blank_database_config_path,
            "submission",
            "populate",
            "--no-confirm",
            s3_metadata.submission_id,
            str(s3_metadata_path),
        ],
    )
    assert result_populate3.exit_code == 0, result_populate3.output
    result_modify3 = runner.invoke(
        cli,
        [
            "db",
            "--config-file",
            blank_database_config_path,
            "submission",
            "modify",
            s3_metadata.submission_id,
            "basic_qc_passed",
            "yes",
        ],
    )
    assert result_modify3.exit_code == 0, result_modify3.output
    result_change1 = runner.invoke(
        cli,
        [
            "db",
            "--config-file",
            blank_database_config_path,
            "submission",
            "change-request",
            s1_metadata.submission_id,
            "delete",
        ],
    )
    assert result_change1.exit_code == 0, result_change1.output

    # generate and check quarterly report
    with runner.isolated_filesystem(temp_dir=tmp_path) as report_tmp_dir:
        result_report = runner.invoke(
            cli,
            ["report", "--config-file", blank_database_config_path, "quarterly", "--year", "2025", "--quarter", "3"],
            catch_exceptions=False,
        )
        assert result_report.exit_code == 0, result_report.output

    overview_output_path = Path(report_tmp_dir) / "Gesamtübersicht.tsv"
    assert overview_output_path.exists()
    with open(overview_output_path, newline="", encoding="utf-8") as overview_file:
        overview_reader = csv.reader(overview_file, delimiter="\t")
        # ignore the header, sort by submitter ID
        rows = sorted(list(overview_reader)[1:], key=itemgetter(3))
    assert rows[0] == [
        "GRZK00007",
        "3",
        "2025",
        "260914050",
        "0",
        "0",
        "2",
        "0",
        "2",
        "0",
        "0",
        "0",
        "1",
        "0",
        "0",
        "1",
    ]
    assert rows[1] == [
        "GRZK00007",
        "3",
        "2025",
        "987654321",
        "1",
        "1",
        "1",
        "1",
        "0",
        "0",
        "1",
        "0",
        "0",
        "0",
        "0",
        "0",
    ]

    dataset_output_path = Path(report_tmp_dir) / "Infos_zu_Datensätzen.tsv"
    assert dataset_output_path.exists()
    with open(dataset_output_path, newline="", encoding="utf-8") as dataset_file:
        dataset_reader = csv.reader(dataset_file, delimiter="\t")
        rows = sorted(list(dataset_reader)[1:], key=itemgetter(3, 4))
    assert rows[0] == [
        "GRZK00007",
        "3",
        "2025",
        "260914050",
        "correction",
        "GKV",
        "oncological",
        "dna",
        "wes",
        "yes",
        "yes",
        "no",
        "NA",
        "not_performed",
        "duo",
        "tumor+germline",
        "index",
        "germline;somatic",
    ]
    assert rows[1] == [
        "GRZK00007",
        "3",
        "2025",
        "260914050",
        "initial",
        "GKV",
        "oncological",
        "dna",
        "wes",
        "yes",
        "yes",
        "yes",
        "NA",
        "not_performed",
        "duo",
        "tumor+germline",
        "index",
        "germline;somatic",
    ]
    assert rows[2] == [
        "GRZK00007",
        "3",
        "2025",
        "987654321",
        "test",
        "GKV",
        "oncological",
        "dna",
        "wes",
        "yes",
        "yes",
        "no",
        "patient refuses to sign consent",
        "no",
        "single",
        "tumor+germline",
        "index",
        "germline;somatic",
    ]

    qc_output_path = Path(report_tmp_dir) / "Detailprüfung.tsv"
    assert qc_output_path.exists()
    with open(qc_output_path, newline="", encoding="utf-8") as qc_file:
        qc_reader = csv.reader(qc_file, delimiter="\t")
        # sort by sequence subtype
        rows = sorted(list(qc_reader)[1:], key=itemgetter(11))
        # one failed submission with two lab datum
        assert len(rows) == 2

    assert rows[0][:13] == [
        "GRZK00007",
        "3",
        "2025",
        "987654321",
        "2025-09-15",
        "test",
        "oncological",
        "single",
        "tumor+germline",
        "index",
        "dna",
        "germline",
        "wes",
    ]


def test_quarterly_migrated_database(blank_database_config_path: Path, tmp_path: Path):
    """Quarterly reports should work on databases migrated from prior schema without backpopulating metadata."""
    # add some minimal test data
    config = DbConfig.from_path(blank_database_config_path)
    tan_g = "a2b6c3d9e8f7123456789abcdef0123456789abcdef0123456789abcdef01234"
    pseudonym = "CASE12345"
    submission_date = "2025-09-14"
    submitter_id = "123456789"
    submission_id = f"{submitter_id}_{submission_date}_d0f805c5"
    with sqlite3.connect(config.db.database_url[len("sqlite:///") :]) as connection:
        connection.execute(
            """
            INSERT INTO submissions(tan_g, pseudonym, id, submission_date, submission_type, submitter_id, data_node_id)
            VALUES(:tan_g, :pseudonym, :id, :submission_date, :submission_type, :submitter_id, :data_node_id)
            """,
            {
                "tan_g": tan_g,
                "pseudonym": pseudonym,
                "id": submission_id,
                "submission_date": submission_date,
                "submission_type": "initial",
                "submitter_id": submitter_id,
                "data_node_id": "GRZXYZ123",
            },
        )

    env = {"GRZ_DB__AUTHOR__PRIVATE_KEY_PASSPHRASE": "test"}

    runner = CliRunner(env=env)
    cli = grzctl.cli.build_cli()

    with runner.isolated_filesystem(temp_dir=tmp_path) as report_tmp_dir:
        result_report = runner.invoke(
            cli,
            ["report", "--config-file", blank_database_config_path, "quarterly", "--year", "2025", "--quarter", "3"],
        )
        assert result_report.exit_code == 0, result_report.output

    overview_output_path = Path(report_tmp_dir) / "Gesamtübersicht.tsv"
    assert overview_output_path.exists()
    with open(overview_output_path, newline="", encoding="utf-8") as overview_file:
        overview_reader = csv.reader(overview_file, delimiter="\t")
        # header + single submitter
        assert len(list(overview_reader)) == 2

    dataset_output_path = Path(report_tmp_dir) / "Infos_zu_Datensätzen.tsv"
    assert dataset_output_path.exists()
    with open(dataset_output_path, newline="", encoding="utf-8") as dataset_file:
        dataset_reader = csv.reader(dataset_file, delimiter="\t")
        # header + single submission
        assert len(list(dataset_reader)) == 2

    qc_output_path = Path(report_tmp_dir) / "Detailprüfung.tsv"
    assert qc_output_path.exists()
    with open(qc_output_path, newline="", encoding="utf-8") as qc_file:
        qc_reader = csv.reader(qc_file, delimiter="\t")
        # header + no detailed QC failures
        assert len(list(qc_reader)) == 1
