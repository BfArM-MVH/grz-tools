"""
Tests for the Prüfbericht submission functionality.
"""

import importlib.resources
import json

import boto3
import click.testing
import grzctl.cli
import pytest
import yaml
from grz_common.workers.upload import S3BotoUploadWorker
from grz_db.models.submission import SubmissionDb, SubmissionStateEnum
from grz_pydantic_models.submission.metadata import (
    REDACTED_LOCAL_CASE_ID,
    REDACTED_TAN,
    GrzSubmissionMetadata,
)

from .. import mock_files
from .common import copy_submission

SUBMISSION_ID = "260914050_2024-07-15_c64603a7"


@pytest.fixture
def archive_db_config_path(tmp_path, migrated_db_config_model):
    """A config with the archives and a migrated DB, whose author key needs no passphrase prompt."""
    data = migrated_db_config_model.model_dump(mode="json", exclude_none=True)
    data["db"]["author"]["private_key_passphrase"] = "test"
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(data))
    return config_path


def test_archive(temp_grzctl_s3_config_file_path, remote_bucket_with_version, working_dir_path, tmp_path):
    submission_dir_ptr = importlib.resources.files(mock_files).joinpath("submissions", "valid_submission")
    with importlib.resources.as_file(submission_dir_ptr) as submission_dir:
        copy_submission(working_dir_path, "encrypted_files", "metadata", source=submission_dir)

        with open(working_dir_path / "metadata" / "metadata.json", mode="r+") as metadata_file:
            metadata_json = json.load(metadata_file)

            # reset donorPseudonym to tanG if index
            for donor in metadata_json["donors"]:
                if donor["relation"] == "index":
                    donor["donorPseudonym"] = metadata_json["submission"]["tanG"]

            # overwrite metadata file
            metadata_file.seek(0)
            json.dump(metadata_json, metadata_file)
            metadata_file.truncate()

        args = [
            "--config",
            temp_grzctl_s3_config_file_path,
            "archive",
            "--submission-dir",
            str(working_dir_path),
            "--no-update-db",
        ]

        runner = click.testing.CliRunner()
        cli = grzctl.cli.build_cli()
        result = runner.invoke(cli, args, catch_exceptions=False)

    # The archive uploads to the "consented" bucket
    consented_bucket = boto3.resource("s3").Bucket("consented")
    uploaded_keys = {o.key for o in consented_bucket.objects.all()}
    assert "260914050_2024-07-15_c64603a7/metadata/metadata.json" in uploaded_keys
    assert "260914050_2024-07-15_c64603a7/logs/progress_upload.cjson" in uploaded_keys
    assert "260914050_2024-07-15_c64603a7/files/target_regions.bed.c4gh" in uploaded_keys

    consented_bucket.download_file(
        Key="260914050_2024-07-15_c64603a7/metadata/metadata.json", Filename=tmp_path / "metadata.json"
    )
    with open(tmp_path / "metadata.json") as metadata_file:
        metadata = GrzSubmissionMetadata.model_validate_json(metadata_file.read())

        # ensure tanG is redacted
        assert metadata.submission.tan_g == REDACTED_TAN

        # ensure local case ID is redacted
        assert metadata.submission.local_case_id == REDACTED_LOCAL_CASE_ID

        # ensure index patient donor pseudonym is redacted
        assert metadata.index_donor.donor_pseudonym == "index"

    assert result.exit_code == 0, result.output


def test_archive_rerun_records_archived(
    archive_db_config_path, migrated_db_config_model, remote_bucket_with_version, working_dir_path, mocker
):
    """A rerun for an archived submission uploads nothing and records ``ARCHIVED``, not an error."""
    submission_dir_ptr = importlib.resources.files(mock_files).joinpath("submissions", "valid_submission")
    with importlib.resources.as_file(submission_dir_ptr) as submission_dir:
        copy_submission(working_dir_path, "encrypted_files", "metadata", source=submission_dir)

    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()
    config_args = ["--config", str(archive_db_config_path)]
    result = runner.invoke(cli, [*config_args, "db", "submission", "add", SUBMISSION_ID], catch_exceptions=False)
    assert result.exit_code == 0, result.output

    archive_args = [*config_args, "archive", "--submission-dir", str(working_dir_path), "--update-db"]
    result = runner.invoke(cli, archive_args, catch_exceptions=False)
    assert result.exit_code == 0, result.output

    upload_spy = mocker.spy(S3BotoUploadWorker, "upload_file")
    result = runner.invoke(cli, archive_args, catch_exceptions=False)
    assert result.exit_code == 0, result.output
    assert upload_spy.call_count == 0, "a rerun must not upload anything"

    db = SubmissionDb(db_url=migrated_db_config_model.db.database_url, author=None)
    states = [state_log.state for state_log in sorted(db.get_submission(SUBMISSION_ID).states, key=lambda s: s.id)]
    assert states == [
        SubmissionStateEnum.ARCHIVING,
        SubmissionStateEnum.ARCHIVED,
        SubmissionStateEnum.ARCHIVING,
        SubmissionStateEnum.ARCHIVED,
    ]
