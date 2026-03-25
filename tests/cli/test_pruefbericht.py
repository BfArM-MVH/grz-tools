"""
Tests for the Prüfbericht submission functionality.
"""

import importlib.resources
import json
import shutil

import click.testing
import grzctl.cli
import pytest
import responses
from grz_pydantic_models.pruefbericht.v0 import LibraryType
from grz_pydantic_models.submission.metadata import REDACTED_TAN

from .. import mock_files

TEST_SUBMISSION_ID = "123456789_1970-01-01_00000000"


@pytest.fixture
def requests_mock(assert_all_requests_are_fired: bool = False):
    with responses.RequestsMock(assert_all_requests_are_fired=assert_all_requests_are_fired) as rsps:
        yield rsps


@pytest.fixture
def bfarm_auth_api(requests_mock):
    """Fakes the endpoint responsible for granting temporary access tokens."""
    requests_mock.post(
        "https://bfarm.localhost/token",
        match=[
            responses.matchers.header_matcher({"Content-Type": "application/x-www-form-urlencoded"}),
            responses.matchers.urlencoded_params_matcher(
                {"grant_type": "client_credentials", "client_id": "pytest", "client_secret": "pysecret"}
            ),
        ],
        json={
            "access_token": "my_token",
            "expires_in": 300,
            "refresh_expires_in": 0,
            "token_type": "Bearer",
            "not-before-policy": 0,
            "scope": "profile email",
        },
    )
    yield requests_mock


@pytest.fixture
def bfarm_submit_api(requests_mock):
    """Fakes the Prüfbericht submission endpoint."""
    # valid submission + valid token
    requests_mock.post(
        "https://bfarm.localhost/api/upload",
        match=[
            responses.matchers.header_matcher({"Authorization": "bearer my_token"}),
            responses.matchers.json_params_matcher(
                {
                    "SubmittedCase": {
                        "submissionDate": "2024-07-15",
                        "submissionType": "test",
                        "tan": "aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000",
                        "submitterId": "260914050",
                        "dataNodeId": "GRZK00007",
                        "diseaseType": "oncological",
                        "dataCategory": "genomic",
                        "libraryType": "wes",
                        "coverageType": "GKV",
                        "dataQualityCheckPassed": True,
                    }
                }
            ),
        ],
    )

    # valid submission + expired token
    requests_mock.post(
        "https://bfarm.localhost/api/upload",
        match=[
            responses.matchers.header_matcher({"Authorization": "bearer expired_token"}),
            responses.matchers.json_params_matcher(
                {
                    "SubmittedCase": {
                        "submissionDate": "2024-07-15",
                        "submissionType": "test",
                        "tan": "aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000",
                        "submitterId": "260914050",
                        "dataNodeId": "GRZK00007",
                        "diseaseType": "oncological",
                        "dataCategory": "genomic",
                        "libraryType": "wes",
                        "coverageType": "GKV",
                        "dataQualityCheckPassed": True,
                    }
                }
            ),
        ],
        status=401,
    )
    yield requests_mock


def test_valid_submission(bfarm_auth_api, bfarm_submit_api, temp_pruefbericht_config_file_path, tmp_path):
    submission_dir_ptr = importlib.resources.files(mock_files).joinpath("submissions", "valid_submission")
    with importlib.resources.as_file(submission_dir_ptr) as submission_dir:
        runner = click.testing.CliRunner(
            env={
                "GRZ_PRUEFBERICHT__AUTHORIZATION_URL": "https://bfarm.localhost/token",
                "GRZ_PRUEFBERICHT__CLIENT_ID": "pytest",
                "GRZ_PRUEFBERICHT__CLIENT_SECRET": "pysecret",
                "GRZ_PRUEFBERICHT__API_BASE_URL": "https://bfarm.localhost/api",
            }
        )
        cli = grzctl.cli.build_cli()

        # generate Prüfbericht JSON
        pruefbericht_json_path = tmp_path / "pruefbericht.json"
        generate_args = ["pruefbericht", "generate", "from-submission-dir", str(submission_dir)]
        generate_result = runner.invoke(cli, generate_args, catch_exceptions=False)
        assert generate_result.exit_code == 0, generate_result.output
        pruefbericht_json_path.write_text(generate_result.output)

        # submit generated Prüfbericht
        submit_args = [
            "pruefbericht",
            "submit",
            "--submission-id",
            TEST_SUBMISSION_ID,
            "--config-file",
            temp_pruefbericht_config_file_path,
            "--pruefbericht-file",
            str(pruefbericht_json_path),
        ]
        submit_result = runner.invoke(cli, submit_args, catch_exceptions=False)

    assert submit_result.exit_code == 0, submit_result.output


def test_valid_submission_with_token(bfarm_submit_api, temp_pruefbericht_config_file_path, tmp_path):
    submission_dir_ptr = importlib.resources.files(mock_files).joinpath("submissions", "valid_submission")
    with importlib.resources.as_file(submission_dir_ptr) as submission_dir:
        runner = click.testing.CliRunner(
            env={
                "GRZ_PRUEFBERICHT__AUTHORIZATION_URL": "https://bfarm.localhost/token",
                "GRZ_PRUEFBERICHT__CLIENT_ID": "pytest",
                "GRZ_PRUEFBERICHT__CLIENT_SECRET": "pysecret",
                "GRZ_PRUEFBERICHT__API_BASE_URL": "https://bfarm.localhost/api",
            }
        )
        cli = grzctl.cli.build_cli()

        # generate Prüfbericht JSON
        pruefbericht_json_path = tmp_path / "pruefbericht.json"
        generate_args = ["pruefbericht", "generate", "from-submission-dir", str(submission_dir)]
        generate_result = runner.invoke(cli, generate_args, catch_exceptions=False)
        assert generate_result.exit_code == 0, generate_result.output
        pruefbericht_json_path.write_text(generate_result.output)

        # submit generated Prüfbericht with a pre-provided token
        submit_args = [
            "pruefbericht",
            "submit",
            "--submission-id",
            TEST_SUBMISSION_ID,
            "--config-file",
            temp_pruefbericht_config_file_path,
            "--pruefbericht-file",
            str(pruefbericht_json_path),
            "--token",
            "my_token",
        ]
        submit_result = runner.invoke(cli, submit_args, catch_exceptions=False)

    assert submit_result.exit_code == 0, submit_result.output


def test_valid_submission_with_expired_token(
    bfarm_auth_api, bfarm_submit_api, temp_pruefbericht_config_file_path, tmp_path
):
    submission_dir_ptr = importlib.resources.files(mock_files).joinpath("submissions", "valid_submission")
    with importlib.resources.as_file(submission_dir_ptr) as submission_dir:
        runner = click.testing.CliRunner(
            env={
                "GRZ_PRUEFBERICHT__AUTHORIZATION_URL": "https://bfarm.localhost/token",
                "GRZ_PRUEFBERICHT__CLIENT_ID": "pytest",
                "GRZ_PRUEFBERICHT__CLIENT_SECRET": "pysecret",
                "GRZ_PRUEFBERICHT__API_BASE_URL": "https://bfarm.localhost/api",
            }
        )
        cli = grzctl.cli.build_cli()

        # generate Prüfbericht JSON
        pruefbericht_json_path = tmp_path / "pruefbericht.json"
        generate_args = ["pruefbericht", "generate", "from-submission-dir", str(submission_dir)]
        generate_result = runner.invoke(cli, generate_args, catch_exceptions=False)
        assert generate_result.exit_code == 0, generate_result.output
        pruefbericht_json_path.write_text(generate_result.output)

        # (try to) submit generated Prüfbericht with an expired token
        submit_args = [
            "pruefbericht",
            "submit",
            "--submission-id",
            TEST_SUBMISSION_ID,
            "--config-file",
            temp_pruefbericht_config_file_path,
            "--pruefbericht-file",
            str(pruefbericht_json_path),
            "--token",
            "expired_token",
        ]
        submit_result = runner.invoke(cli, submit_args, catch_exceptions=False)

    assert submit_result.exit_code == 0, submit_result.output


def test_generate_pruefbericht_multiple_library_types(temp_pruefbericht_config_file_path, tmp_path):
    submission_dir_ptr = importlib.resources.files(mock_files).joinpath("submissions", "valid_submission")
    with importlib.resources.as_file(submission_dir_ptr) as submission_dir:
        # create and modify a temporary copy of the metadata JSON
        shutil.copytree(submission_dir, tmp_path, dirs_exist_ok=True)
        with open(tmp_path / "metadata" / "metadata.json", mode="r+") as metadata_file:
            metadata = json.load(metadata_file)

            # use to differentiate from standard valid submission mock
            metadata["submission"]["tanG"] = "aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000001"
            # set only ONE libary type to WGS
            # should sort higher than the other WES data and therefore WGS is sent in Pruefbericht
            metadata["donors"][0]["labData"][0]["libraryType"] = "wgs"
            metadata["donors"][0]["labData"][0]["sequenceData"]["minCoverage"] = 20

            metadata_file.seek(0)
            json.dump(metadata, metadata_file)
            metadata_file.truncate()

        args = [
            "pruefbericht",
            "generate",
            "from-submission-dir",
            str(tmp_path),
        ]

        runner = click.testing.CliRunner()
        cli = grzctl.cli.build_cli()
        result = runner.invoke(cli, args, catch_exceptions=False)

    assert result.exit_code == 0, result.output
    # Check that the generated Pruefbericht correctly selected the most expensive library type
    pruefbericht_data = json.loads(result.output)
    assert pruefbericht_data["SubmittedCase"]["libraryType"] == "wgs"


def test_generate_fails_with_invalid_library_type(temp_pruefbericht_config_file_path, tmp_path):
    submission_dir_ptr = importlib.resources.files(mock_files).joinpath("submissions", "valid_submission")
    with importlib.resources.as_file(submission_dir_ptr) as submission_dir:
        # create and modify a temporary copy of the metadata JSON
        shutil.copytree(submission_dir, tmp_path, dirs_exist_ok=True)
        with open(tmp_path / "metadata" / "metadata.json", mode="r+") as metadata_file:
            metadata = json.load(metadata_file)

            # remove other donors/lab data
            metadata["donors"] = metadata["donors"][:1]
            metadata["donors"][0]["labData"] = metadata["donors"][0]["labData"][:1]
            metadata["submission"]["genomicStudyType"] = "single"
            # set to valid submission library type but invalid pruefbericht library type
            metadata["donors"][0]["labData"][0]["libraryType"] = "other"

            metadata_file.seek(0)
            json.dump(metadata, metadata_file)
            metadata_file.truncate()

        args = [
            "pruefbericht",
            "generate",
            "from-submission-dir",
            str(tmp_path),
        ]

        runner = click.testing.CliRunner()
        cli = grzctl.cli.build_cli()
        result = runner.invoke(cli, args, catch_exceptions=False)
        assert result.exit_code != 0, result.output


def test_refuse_redacted_tang(temp_pruefbericht_config_file_path, tmp_path):
    submission_dir_ptr = importlib.resources.files(mock_files).joinpath("submissions", "valid_submission")
    with importlib.resources.as_file(submission_dir_ptr) as submission_dir:
        # create and modify a temporary copy of the metadata JSON
        shutil.copytree(submission_dir, tmp_path, dirs_exist_ok=True)
        with open(tmp_path / "metadata" / "metadata.json", mode="r+") as metadata_file:
            metadata = json.load(metadata_file)

            # set tanG to REDACTED_TAN
            metadata["submission"]["tanG"] = REDACTED_TAN

            metadata_file.seek(0)
            json.dump(metadata, metadata_file)
            metadata_file.truncate()

        runner = click.testing.CliRunner(
            env={
                "GRZ_PRUEFBERICHT__AUTHORIZATION_URL": "https://bfarm.localhost/token",
                "GRZ_PRUEFBERICHT__CLIENT_ID": "pytest",
                "GRZ_PRUEFBERICHT__CLIENT_SECRET": "pysecret",
                "GRZ_PRUEFBERICHT__API_BASE_URL": "https://bfarm.localhost/api",
            }
        )
        cli = grzctl.cli.build_cli()

        # generate Prüfbericht JSON
        pruefbericht_json_path = tmp_path / "pruefbericht.json"
        generate_args = ["pruefbericht", "generate", "from-submission-dir", str(tmp_path)]
        generate_result = runner.invoke(cli, generate_args, catch_exceptions=False)
        assert generate_result.exit_code == 0, generate_result.output
        pruefbericht_json_path.write_text(generate_result.output)

        # attempt to submit
        submit_args = [
            "pruefbericht",
            "submit",
            "--submission-id",
            TEST_SUBMISSION_ID,
            "--config-file",
            temp_pruefbericht_config_file_path,
            "--pruefbericht-file",
            str(pruefbericht_json_path),
        ]
        with pytest.raises(ValueError, match=r"Refusing to submit a Prüfbericht with a redacted TAN"):
            runner.invoke(cli, submit_args, catch_exceptions=False)


@pytest.fixture
def pruefbericht_db_config(tmp_path):
    """Create a test database config for pruefbericht tests."""
    import json

    db_path = tmp_path / "test.db"
    db_url = f"sqlite:///{db_path}"

    config = {"db": {"database_url": db_url, "author": {"name": "test_author"}}}

    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config))

    return {"db_path": db_path, "db_url": db_url, "config": config, "config_path": config_path}


def test_generate_from_database(temp_pruefbericht_config_file_path, pruefbericht_db_config):
    """Test generating Prüfbericht from database."""
    import datetime
    import json

    from grz_db.models.submission import Donor, SubmissionDb
    from grz_pydantic_models.submission.metadata import (
        LibraryType,
        Relation,
        SequenceSubtype,
        SequenceType,
    )
    from grz_pydantic_models.submission.metadata.v1 import (
        CoverageType,
        DiseaseType,
        SubmissionType,
    )

    # Use the fixture
    db_url = pruefbericht_db_config["db_url"]
    config_path = pruefbericht_db_config["config_path"]

    db = SubmissionDb(db_url=db_url, author=None, debug=False)
    db.initialize_schema()

    db.add_submission(TEST_SUBMISSION_ID)
    # Populate submission fields to match the mock data used in other tests
    db.modify_submission(TEST_SUBMISSION_ID, "submission_date", datetime.date(2024, 7, 15))
    db.modify_submission(TEST_SUBMISSION_ID, "submission_type", SubmissionType.test)
    db.modify_submission(
        TEST_SUBMISSION_ID, "tan_g", "aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000"
    )
    db.modify_submission(TEST_SUBMISSION_ID, "submitter_id", "260914050")
    db.modify_submission(TEST_SUBMISSION_ID, "data_node_id", "GRZK00007")
    db.modify_submission(TEST_SUBMISSION_ID, "disease_type", DiseaseType.oncological)
    db.modify_submission(TEST_SUBMISSION_ID, "coverage_type", CoverageType.GKV)

    # Add index donor with library types
    index_donor = Donor(
        submission_id=TEST_SUBMISSION_ID,
        pseudonym="test_donor",
        relation=Relation.index_,
        library_types={LibraryType.wes, LibraryType.panel},
        sequence_types={SequenceType.dna},
        sequence_subtypes={SequenceSubtype.germline},
        mv_consented=True,
        research_consented=True,
    )
    db.add_donor(index_donor)

    # Test the from-database command
    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()

    args = [
        "pruefbericht",
        "generate",
        "from-database",
        "--submission-id",
        TEST_SUBMISSION_ID,
        "--config-file",
        str(config_path),
    ]

    result = runner.invoke(cli, args, catch_exceptions=False)

    assert result.exit_code == 0, result.output

    # Verify the generated Prüfbericht matches expected values
    pruefbericht_data = json.loads(result.output)
    assert pruefbericht_data["SubmittedCase"] == {
        "submissionDate": "2024-07-15",
        "submissionType": "test",
        "tan": "aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000",
        "submitterId": "260914050",
        "dataNodeId": "GRZK00007",
        "diseaseType": "oncological",
        "dataCategory": "genomic",
        "libraryType": "wes",  # Most expensive
        "coverageType": "GKV",
        "dataQualityCheckPassed": True,
    }


def test_generate_from_database_missing_submission(pruefbericht_db_config):
    """Test error when submission doesn't exist in database."""
    from grz_db.models.submission import SubmissionDb

    # Use the fixture
    db_url = pruefbericht_db_config["db_url"]
    config_path = pruefbericht_db_config["config_path"]

    db = SubmissionDb(db_url=db_url, author=None, debug=False)
    db.initialize_schema()

    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()

    args = [
        "pruefbericht",
        "generate",
        "from-database",
        "--submission-id",
        "nonexistent_id",
        "--config-file",
        str(config_path),
    ]

    result = runner.invoke(cli, args, catch_exceptions=False)

    assert result.exit_code != 0
    assert "not found in database" in result.output


def test_generate_from_database_missing_fields(pruefbericht_db_config):
    """Test error when submission has missing required fields."""
    from grz_db.models.submission import SubmissionDb

    # Use the fixture
    db_url = pruefbericht_db_config["db_url"]
    config_path = pruefbericht_db_config["config_path"]

    db = SubmissionDb(db_url=db_url, author=None, debug=False)
    db.initialize_schema()

    # Add submission without populating required fields
    db.add_submission(TEST_SUBMISSION_ID)

    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()

    args = [
        "pruefbericht",
        "generate",
        "from-database",
        "--submission-id",
        TEST_SUBMISSION_ID,
        "--config-file",
        str(config_path),
    ]

    result = runner.invoke(cli, args, catch_exceptions=False)

    assert result.exit_code != 0
    assert "missing required fields" in result.output


def test_generate_from_database_no_index_donor(pruefbericht_db_config):
    """Test error when no index donor exists."""
    import datetime

    from grz_db.models.submission import Donor, SubmissionDb
    from grz_pydantic_models.submission.metadata import (
        LibraryType,
        Relation,
        SequenceSubtype,
        SequenceType,
    )
    from grz_pydantic_models.submission.metadata.v1 import (
        CoverageType,
        DiseaseType,
        SubmissionType,
    )

    # Use the fixture
    db_url = pruefbericht_db_config["db_url"]
    config_path = pruefbericht_db_config["config_path"]

    db = SubmissionDb(db_url=db_url, author=None, debug=False)
    db.initialize_schema()
    db.add_submission(TEST_SUBMISSION_ID)

    # Populate all required fields (matching the pattern from test_generate_from_database)
    db.modify_submission(TEST_SUBMISSION_ID, "submission_date", datetime.date(2024, 7, 15))
    db.modify_submission(TEST_SUBMISSION_ID, "submission_type", SubmissionType.test)
    db.modify_submission(
        TEST_SUBMISSION_ID, "tan_g", "aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000"
    )
    db.modify_submission(TEST_SUBMISSION_ID, "submitter_id", "260914050")
    db.modify_submission(TEST_SUBMISSION_ID, "data_node_id", "GRZK00007")
    db.modify_submission(TEST_SUBMISSION_ID, "disease_type", DiseaseType.oncological)
    db.modify_submission(TEST_SUBMISSION_ID, "coverage_type", CoverageType.GKV)

    # Add a non-index donor such as mother
    donor = Donor(
        submission_id=TEST_SUBMISSION_ID,
        pseudonym="test_donor",
        relation=Relation.mother,  # NOT index
        library_types={LibraryType.wes},
        sequence_types={SequenceType.dna},
        sequence_subtypes={SequenceSubtype.germline},
        mv_consented=True,
        research_consented=True,
    )
    db.add_donor(donor)

    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()

    args = [
        "pruefbericht",
        "generate",
        "from-database",
        "--submission-id",
        TEST_SUBMISSION_ID,
        "--config-file",
        str(config_path),
    ]

    result = runner.invoke(cli, args, catch_exceptions=False)

    assert result.exit_code != 0
    assert "No index donor found" in result.output


class TestLibraryTypeMostExpensive:
    def test_single_type_returned(self):
        assert LibraryType.most_expensive({"wes"}) == LibraryType.wes

    def test_most_expensive_is_wgs_lr(self):
        assert LibraryType.most_expensive({"panel", "wes", "wgs", "wgs_lr"}) == LibraryType.wgs_lr

    def test_most_expensive_excludes_none(self):
        """'none' must never win over a real sequencing type."""
        assert LibraryType.most_expensive({"panel", "none"}) == LibraryType.panel

    def test_none_alone_is_valid(self):
        assert LibraryType.most_expensive({"none"}) == LibraryType.none

    def test_invalid_types_are_ignored(self):
        assert LibraryType.most_expensive({"wgs", "invalid_type"}) == LibraryType.wgs

    def test_all_invalid_raises(self):
        with pytest.raises(ValueError, match="cannot be submitted"):
            LibraryType.most_expensive({"invalid_type", "another_bad_one"})

    def test_priority_ordering(self):
        """Verify the full priority ladder explicitly."""
        ordered = sorted(LibraryType, key=lambda lt: lt.reimbursement_priority)
        assert ordered == [
            LibraryType.none,
            LibraryType.panel,
            LibraryType.wes,
            LibraryType.wgs,
            LibraryType.wgs_lr,
        ]
