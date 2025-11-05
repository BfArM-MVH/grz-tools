"""Tests for correctly failing metadata."""

import importlib.resources

import pytest
from grz_pydantic_models.submission.metadata.v1 import GrzSubmissionMetadata
from pydantic import ValidationError

from . import resources

resource_files = importlib.resources.files(resources)

metadata_missing_read_order = resource_files.joinpath("failing_metadata/missing-read-order.json")
metadata_missing_vcf_file = resource_files.joinpath("failing_metadata/missing-vcf-file.json")
metadata_missing_fastq_r2 = resource_files.joinpath("failing_metadata/missing-fastq-r2.json")
metadata_no_target_regions = resource_files.joinpath("failing_metadata/missing-target-regions.json")
metadata_incompatible_reference_genomes = resource_files.joinpath(
    "failing_metadata/incompatible-reference-genomes.json"
)
metadata_duplicate_run_id = resource_files.joinpath("failing_metadata/duplicate-run-id.json")


def test_submission_metadata_fails():
    error_types = (ValueError, ValidationError, SystemExit)
    with pytest.raises(error_types, match="No read order specified for FASTQ file"):
        GrzSubmissionMetadata.model_validate_json(metadata_missing_read_order.read_text())

    with pytest.raises(error_types, match="BED file missing for lab datum"):
        GrzSubmissionMetadata.model_validate_json(metadata_no_target_regions.read_text())

    # missing VCF is allowed
    GrzSubmissionMetadata.model_validate_json(metadata_missing_vcf_file.read_text())

    with pytest.raises(
        error_types, match="Paired end sequencing layout but not there is not exactly one R1 and one R2"
    ):
        GrzSubmissionMetadata.model_validate_json(metadata_missing_fastq_r2.read_text())

    with pytest.raises(error_types, match="Incompatible reference genomes found"):
        GrzSubmissionMetadata.model_validate_json(metadata_incompatible_reference_genomes.read_text())

    with pytest.raises(error_types, match="must have a unique combination of flowcell_id, lane_id, and read_order"):
        GrzSubmissionMetadata.model_validate_json(metadata_duplicate_run_id.read_text())
