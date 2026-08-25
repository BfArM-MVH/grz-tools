"""
Tests for the Prüfbericht submission functionality.
"""

import importlib.resources
import json

import boto3
import click.testing
import grzctl.cli
from grz_pydantic_models.submission.metadata import (
    REDACTED_LOCAL_CASE_ID,
    REDACTED_TAN,
    GrzSubmissionMetadata,
)

from .. import mock_files
from .common import copy_submission


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
