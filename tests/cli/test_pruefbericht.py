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
from grz_pydantic_models.submission.metadata import REDACTED_TAN

from .. import mock_files


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
            "--config-file",
            temp_pruefbericht_config_file_path,
            "--pruefbericht-file",
            str(pruefbericht_json_path),
        ]
        with pytest.raises(ValueError, match="Refusing to submit a Prüfbericht with a redacted TAN"):
            runner.invoke(cli, submit_args, catch_exceptions=False)
