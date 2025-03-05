"""Tests for the parser module."""

from pathlib import Path

import pytest
from pydantic import ValidationError

from grz_cli.parser import EncryptedSubmission, SubmissionMetadata

metadata_missing_read_order = "tests/mock_files/metadata_validation/missing-read-order.json"
metadata_missing_vcf_file = "tests/mock_files/metadata_validation/missing-vcf-file.json"
metadata_missing_fastq_r2 = "tests/mock_files/metadata_validation/missing-fastq-r2.json"
metadata_no_target_regions = "tests/mock_files/metadata_validation/missing-target-regions.json"
metadata_incompatible_reference_genomes = "tests/mock_files/metadata_validation/incompatible-reference-genomes.json"


def test_submission_metadata(temp_metadata_file_path):
    submission_metadata = SubmissionMetadata(temp_metadata_file_path)

    errors = list(submission_metadata.validate())
    assert errors == []

    assert len(submission_metadata.files) > 0


def test_submission_metadata_fails():
    error_types = (ValueError, ValidationError, SystemExit)
    with pytest.raises(error_types, match="No read order specified for FASTQ file"):
        SubmissionMetadata(metadata_missing_read_order)

    with pytest.raises(error_types, match="BED file missing for lab datum"):
        SubmissionMetadata(metadata_no_target_regions)

    with pytest.raises(error_types, match="VCF file missing for lab datum"):
        SubmissionMetadata(metadata_missing_vcf_file)

    with pytest.raises(error_types, match="Paired end sequencing layout but missing R2 file for flowcell id"):
        SubmissionMetadata(metadata_missing_fastq_r2)

    with pytest.raises(error_types, match="Incompatible reference genomes found"):
        SubmissionMetadata(metadata_incompatible_reference_genomes)


def test_encrypted_submission():
    input_path = "/submission/files/a.fastq"

    for i in (input_path, Path(input_path)):
        enc_path = EncryptedSubmission.get_encrypted_file_path(i)
        assert enc_path == Path(input_path + ".c4gh")
        enc_path = EncryptedSubmission.get_encryption_header_path(input_path)
        assert enc_path == Path(input_path + ".c4gh_header")


# TODO: test encrypt submission
# TODO: test upload submission
