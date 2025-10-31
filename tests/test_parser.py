"""Tests for the parser module."""

from pathlib import Path

from grz_common.workers.submission import EncryptedSubmission, SubmissionMetadata


def test_submission_metadata(temp_metadata_file_path, identifiers_config_model):
    submission_metadata = SubmissionMetadata(temp_metadata_file_path)

    errors = list(submission_metadata.validate(identifiers_config_model.identifiers))
    assert errors == []

    assert len(submission_metadata.files) > 0


def test_encrypted_submission():
    input_path = "/submission/files/a.fastq"

    for i in (input_path, Path(input_path)):
        enc_path = EncryptedSubmission.get_encrypted_file_path(i)
        assert enc_path == Path(input_path + ".c4gh")
        enc_path = EncryptedSubmission.get_encryption_header_path(input_path)
        assert enc_path == Path(input_path + ".c4gh_header")
