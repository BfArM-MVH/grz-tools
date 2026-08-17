"""Tests for correctly failing metadata."""

import importlib.resources

import pytest
from grz_pydantic_models.submission.metadata.v1 import GrzSubmissionMetadata
from grz_pydantic_models_testing import failing_metadata
from pydantic import ValidationError

resource_files = importlib.resources.files(failing_metadata)

error_types = (ValueError, ValidationError, SystemExit)


@pytest.mark.parametrize(
    "fixture_name,expected_match",
    (
        ("missing-read-order.json", "No read order specified for FASTQ file"),
        ("missing-target-regions.json", "BED file missing for lab datum"),
        ("missing-fastq-r2.json", "Paired end sequencing layout but not there is not exactly one R1 and one R2"),
        ("incompatible-reference-genomes.json", "Incompatible reference genomes found"),
        ("duplicate-run-id.json", "must have a unique combination of flowcell_id, lane_id, and read_order"),
    ),
)
def test_submission_metadata_fails(fixture_name: str, expected_match: str):
    with pytest.raises(error_types, match=expected_match):
        GrzSubmissionMetadata.model_validate_json(resource_files.joinpath(fixture_name).read_text())


def test_submission_metadata_missing_vcf_file_passes():
    """A missing VCF file is allowed."""
    GrzSubmissionMetadata.model_validate_json(resource_files.joinpath("missing-vcf-file.json").read_text())
