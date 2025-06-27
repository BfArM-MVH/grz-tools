"""
Tests for the Prüfbericht submission functionality.
"""

import importlib.resources
import json
import shutil

import click.testing
import grzctl
from grz_pydantic_models.submission.metadata import GrzSubmissionMetadata, Relation

from .. import mock_files


def test_archive(temp_s3_config_file_path, remote_bucket, working_dir_path, tmp_path):
    submission_dir_ptr = importlib.resources.files(mock_files).joinpath("submissions", "valid_submission")
    with importlib.resources.as_file(submission_dir_ptr) as submission_dir:
        shutil.copytree(submission_dir / "encrypted_files", working_dir_path / "encrypted_files", dirs_exist_ok=True)
        shutil.copytree(submission_dir / "metadata", working_dir_path / "metadata", dirs_exist_ok=True)

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
            "archive",
            "--config-file",
            temp_s3_config_file_path,
            "--submission-dir",
            str(working_dir_path),
        ]

        runner = click.testing.CliRunner(
            env={
                "GRZ_S3_OPTIONS__ENDPOINT_URL": "",
            }
        )
        cli = grzctl.cli.build_cli()
        result = runner.invoke(cli, args, catch_exceptions=False)

    uploaded_keys = {o.key for o in remote_bucket.objects.all()}
    assert "260914050_2024-07-15_c64603a7/metadata/metadata.json" in uploaded_keys
    assert "260914050_2024-07-15_c64603a7/logs/progress_upload.cjson" in uploaded_keys
    assert "260914050_2024-07-15_c64603a7/files/target_regions.bed.c4gh" in uploaded_keys

    remote_bucket.download_file(
        Key="260914050_2024-07-15_c64603a7/metadata/metadata.json", Filename=tmp_path / "metadata.json"
    )
    with open(tmp_path / "metadata.json") as metadata_file:
        metadata = GrzSubmissionMetadata.model_validate_json(metadata_file.read())

        # ensure tanG is redacted
        assert metadata.submission.tan_g == "".join(["0"] * 64)

        # ensure local case ID is redacted
        assert not metadata.submission.local_case_id

        # ensure index patient donor pseudonym is redacted
        index_patient = next(donor for donor in metadata.donors if donor.relation == Relation.index_)
        assert index_patient.donor_pseudonym == "index"

    assert result.exit_code == 0, result.output
