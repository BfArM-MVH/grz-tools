"""
Tests for the Prüfbericht submission functionality.
"""

import datetime
import importlib.resources
import json
import shutil

import click.testing
import grzctl.cli
import pytest
import responses

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
                        "submissionType": "initial",
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

    # valid submission with multiple library types + valid token
    requests_mock.post(
        "https://bfarm.localhost/api/upload",
        match=[
            responses.matchers.header_matcher({"Authorization": "bearer my_token"}),
            responses.matchers.json_params_matcher(
                {
                    "SubmittedCase": {
                        "submissionDate": "2024-07-15",
                        "submissionType": "initial",
                        "tan": "aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000000aaaaaaaa00000001",
                        "submitterId": "260914050",
                        "dataNodeId": "GRZK00007",
                        "diseaseType": "oncological",
                        "dataCategory": "genomic",
                        "libraryType": "wgs",
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
                        "submissionType": "initial",
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

    # valid test submission + valid token
    requests_mock.post(
        "https://bfarm.localhost/api/upload",
        match=[
            responses.matchers.header_matcher({"Authorization": "bearer my_token"}),
            responses.matchers.json_params_matcher(
                {
                    "SubmittedCase": {
                        "submissionDate": "2001-01-01",
                        "submissionType": "initial",
                        "tan": "TESTGRZXYZ123123456789C05D8A6B890A7B84E9A2F3C1D086FE24C05D8A6B89",
                        "submitterId": "123456789",
                        "dataNodeId": "GRZXYZ123",
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

    yield requests_mock


def test_valid_submission(bfarm_auth_api, bfarm_submit_api, temp_pruefbericht_config_file_path):
    submission_dir_ptr = importlib.resources.files(mock_files).joinpath("submissions", "valid_submission")
    with importlib.resources.as_file(submission_dir_ptr) as submission_dir:
        args = [
            "pruefbericht",
            "--config-file",
            temp_pruefbericht_config_file_path,
            "--submission-dir",
            str(submission_dir),
        ]

        runner = click.testing.CliRunner(
            env={
                "GRZ_PRUEFBERICHT__AUTHORIZATION_URL": "https://bfarm.localhost/token",
                "GRZ_PRUEFBERICHT__CLIENT_ID": "pytest",
                "GRZ_PRUEFBERICHT__CLIENT_SECRET": "pysecret",
                "GRZ_PRUEFBERICHT__API_BASE_URL": "https://bfarm.localhost/api",
            }
        )
        cli = grzctl.cli.build_cli()
        result = runner.invoke(cli, args, catch_exceptions=False)

    assert result.exit_code == 0, result.output


def test_valid_submission_with_json_output(bfarm_auth_api, bfarm_submit_api, temp_pruefbericht_config_file_path):
    submission_dir_ptr = importlib.resources.files(mock_files).joinpath("submissions", "valid_submission")
    with importlib.resources.as_file(submission_dir_ptr) as submission_dir:
        args = [
            "pruefbericht",
            "--config-file",
            temp_pruefbericht_config_file_path,
            "--submission-dir",
            str(submission_dir),
            "--json",
        ]

        runner = click.testing.CliRunner(
            env={
                "GRZ_PRUEFBERICHT__AUTHORIZATION_URL": "https://bfarm.localhost/token",
                "GRZ_PRUEFBERICHT__CLIENT_ID": "pytest",
                "GRZ_PRUEFBERICHT__CLIENT_SECRET": "pysecret",
                "GRZ_PRUEFBERICHT__API_BASE_URL": "https://bfarm.localhost/api",
            }
        )
        cli = grzctl.cli.build_cli()
        result = runner.invoke(cli, args, catch_exceptions=False)

    assert result.exit_code == 0, result.output

    output = json.loads(result.output)
    datetime.datetime.fromisoformat(output["expires"])
    assert output["token"] == "my_token"


def test_valid_submission_with_token(bfarm_submit_api, temp_pruefbericht_config_file_path):
    submission_dir_ptr = importlib.resources.files(mock_files).joinpath("submissions", "valid_submission")
    with importlib.resources.as_file(submission_dir_ptr) as submission_dir:
        args = [
            "pruefbericht",
            "--config-file",
            temp_pruefbericht_config_file_path,
            "--submission-dir",
            str(submission_dir),
            "--token",
            "my_token",
        ]

        runner = click.testing.CliRunner(
            env={
                "GRZ_PRUEFBERICHT__AUTHORIZATION_URL": "https://bfarm.localhost/token",
                "GRZ_PRUEFBERICHT__CLIENT_ID": "pytest",
                "GRZ_PRUEFBERICHT__CLIENT_SECRET": "pysecret",
                "GRZ_PRUEFBERICHT__API_BASE_URL": "https://bfarm.localhost/api",
            }
        )
        cli = grzctl.cli.build_cli()
        result = runner.invoke(cli, args, catch_exceptions=False)

    assert result.exit_code == 0, result.output


def test_valid_submission_with_expired_token(bfarm_auth_api, bfarm_submit_api, temp_pruefbericht_config_file_path):
    submission_dir_ptr = importlib.resources.files(mock_files).joinpath("submissions", "valid_submission")
    with importlib.resources.as_file(submission_dir_ptr) as submission_dir:
        args = [
            "pruefbericht",
            "--config-file",
            temp_pruefbericht_config_file_path,
            "--submission-dir",
            str(submission_dir),
            "--token",
            "expired_token",
        ]

        runner = click.testing.CliRunner(
            env={
                "GRZ_PRUEFBERICHT__AUTHORIZATION_URL": "https://bfarm.localhost/token",
                "GRZ_PRUEFBERICHT__CLIENT_ID": "pytest",
                "GRZ_PRUEFBERICHT__CLIENT_SECRET": "pysecret",
                "GRZ_PRUEFBERICHT__API_BASE_URL": "https://bfarm.localhost/api",
            }
        )
        cli = grzctl.cli.build_cli()
        result = runner.invoke(cli, args, catch_exceptions=False)

    assert result.exit_code == 0, result.output


def test_valid_submission_multiple_library_types(
    bfarm_auth_api, bfarm_submit_api, temp_pruefbericht_config_file_path, tmp_path
):
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

            metadata_file.seek(0)
            json.dump(metadata, metadata_file)
            metadata_file.truncate()

        args = [
            "pruefbericht",
            "--config-file",
            temp_pruefbericht_config_file_path,
            "--submission-dir",
            str(tmp_path),
        ]

        runner = click.testing.CliRunner(
            env={
                "GRZ_PRUEFBERICHT__AUTHORIZATION_URL": "https://bfarm.localhost/token",
                "GRZ_PRUEFBERICHT__CLIENT_ID": "pytest",
                "GRZ_PRUEFBERICHT__CLIENT_SECRET": "pysecret",
                "GRZ_PRUEFBERICHT__API_BASE_URL": "https://bfarm.localhost/api",
            }
        )
        cli = grzctl.cli.build_cli()
        result = runner.invoke(cli, args, catch_exceptions=False)

    assert result.exit_code == 0, result.output


def test_invalid_submission_invalid_library_type(
    bfarm_auth_api, bfarm_submit_api, temp_pruefbericht_config_file_path, tmp_path
):
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
            "--config-file",
            temp_pruefbericht_config_file_path,
            "--submission-dir",
            str(tmp_path),
        ]

        runner = click.testing.CliRunner(
            env={
                "GRZ_PRUEFBERICHT__AUTHORIZATION_URL": "https://bfarm.localhost/token",
                "GRZ_PRUEFBERICHT__CLIENT_ID": "pytest",
                "GRZ_PRUEFBERICHT__CLIENT_SECRET": "pysecret",
                "GRZ_PRUEFBERICHT__API_BASE_URL": "https://bfarm.localhost/api",
            }
        )
        cli = grzctl.cli.build_cli()
        with pytest.raises(ValueError, match="cannot be submitted in the Prüfbericht"):
            runner.invoke(cli, args, catch_exceptions=False)


def test_valid_test_submission(bfarm_auth_api, bfarm_submit_api, temp_pruefbericht_config_file_path, tmp_path):
    submission_dir_ptr = importlib.resources.files(mock_files).joinpath("submissions", "valid_submission")
    with importlib.resources.as_file(submission_dir_ptr) as submission_dir:
        # create and modify a temporary copy of the metadata JSON
        shutil.copytree(submission_dir, tmp_path, dirs_exist_ok=True)
        with open(tmp_path / "metadata" / "metadata.json", mode="r+") as metadata_file:
            metadata = json.load(metadata_file)

            metadata["submission"]["submissionDate"] = "2001-01-01"
            metadata["submission"]["tanG"] = "TESTGRZXYZ123123456789C05D8A6B890A7B84E9A2F3C1D086FE24C05D8A6B89"
            metadata["submission"]["genomicDataCenterId"] = "GRZXYZ123"
            metadata["submission"]["submitterId"] = "123456789"

            metadata_file.seek(0)
            json.dump(metadata, metadata_file)
            metadata_file.truncate()

        args = [
            "pruefbericht",
            "--config-file",
            temp_pruefbericht_config_file_path,
            "--submission-dir",
            str(tmp_path),
        ]

        runner = click.testing.CliRunner(
            env={
                "GRZ_PRUEFBERICHT__AUTHORIZATION_URL": "https://bfarm.localhost/token",
                "GRZ_PRUEFBERICHT__CLIENT_ID": "pytest",
                "GRZ_PRUEFBERICHT__CLIENT_SECRET": "pysecret",
                "GRZ_PRUEFBERICHT__API_BASE_URL": "https://bfarm.localhost/api",
            }
        )
        cli = grzctl.cli.build_cli()

        # submitting Prüfbericht for test submission should fail without --test
        result_fail = runner.invoke(cli, args, catch_exceptions=False)
        assert result_fail.exit_code != 0, result_fail.output

        result_pass = runner.invoke(cli, [*args, "--test"], catch_exceptions=False)
        assert result_pass.exit_code == 0, result_pass.output
