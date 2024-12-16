"""Tests for the parser module."""

from pathlib import Path

import pytest

from grz_cli.parser import EncryptedSubmission, SubmissionMetadata

metadata_missing_read_order = "tests/mock_files/metadata_validation/missing-read-order.json"
metadata_missing_vcf_file = "tests/mock_files/metadata_validation/missing-vcf-file.json"
metadata_missing_fastq_r2 = "tests/mock_files/metadata_validation/missing-fastq-r2.json"
metadata_no_target_regions = "tests/mock_files/metadata_validation/missing-target-regions.json"


def test_submission_metadata(temp_metadata_file_path):
    submission_metadata = SubmissionMetadata(temp_metadata_file_path)

    errors = list(submission_metadata.validate())
    assert errors == []

    assert len(submission_metadata.files) > 0


def test_submission_metadata_fails():
    with pytest.raises(ValueError, match="No read order specified for FASTQ file"):
        SubmissionMetadata(metadata_missing_read_order)

    with pytest.raises(ValueError, match="BED file missing for lab datum"):
        SubmissionMetadata(metadata_no_target_regions)

    with pytest.raises(ValueError, match="VCF file missing for lab datum"):
        SubmissionMetadata(metadata_missing_vcf_file)

    with pytest.raises(ValueError, match="Paired end sequencing layout but number of"):
        SubmissionMetadata(metadata_missing_fastq_r2)


def test_encrypted_submission():
    input_path = "/submission/files/a.fastq"

    for i in (input_path, Path(input_path)):
        enc_path = EncryptedSubmission.get_encrypted_file_path(i)
        assert enc_path == Path(input_path + ".c4gh")
        enc_path = EncryptedSubmission.get_encryption_header_path(input_path)
        assert enc_path == Path(input_path + ".c4gh_header")


# TODO: test encrypt submission
# TODO: test upload submission
