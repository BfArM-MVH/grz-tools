"""Tests for the parser module."""

import json
from pathlib import Path

from grz_common.workers.submission import EncryptedSubmission, SubmissionMetadata


def test_submission_metadata(temp_metadata_file_path, identifiers_config_model):
    submission_metadata = SubmissionMetadata(temp_metadata_file_path)

    errors = list(submission_metadata.validate(identifiers_config_model.identifiers))
    assert errors == []

    assert len(submission_metadata.files) > 0


def test_submission_metadata_reports_a_consent_datetime_without_a_timezone(
    temp_metadata_file_path, identifiers_config_model
):
    """
    FHIR requires a timezone alongside a time, so a new submission may not state one without.

    The models keep such a value so that documents written before the rule was enforced stay
    readable, which means the document still parses. Refusing it is a submission-time rule, so it
    surfaces as a validation error naming the offending value rather than as a parse failure.
    """
    metadata = json.loads(temp_metadata_file_path.read_text())
    scope = metadata["donors"][0]["researchConsents"][0]["scope"]
    assert scope["provision"]["period"]["start"], "fixture must carry a consent period to amend"
    scope["provision"]["period"]["start"] = "2020-09-01T14:37:22"
    temp_metadata_file_path.write_text(json.dumps(metadata))

    errors = list(SubmissionMetadata(temp_metadata_file_path).validate(identifiers_config_model.identifiers))

    assert [error for error in errors if "2020-09-01T14:37:22" in error and "timezone" in error], errors


def test_submission_metadata_accepts_a_consent_datetime_with_a_timezone(
    temp_metadata_file_path, identifiers_config_model
):
    """The same value is accepted once it states a zone, so the check is about the zone alone."""
    metadata = json.loads(temp_metadata_file_path.read_text())
    metadata["donors"][0]["researchConsents"][0]["scope"]["provision"]["period"]["start"] = "2020-09-01T14:37:22+02:00"
    temp_metadata_file_path.write_text(json.dumps(metadata))

    submission_metadata = SubmissionMetadata(temp_metadata_file_path)

    assert list(submission_metadata.validate(identifiers_config_model.identifiers)) == []


def test_encrypted_submission():
    input_path = "/submission/files/a.fastq"

    for i in (input_path, Path(input_path)):
        enc_path = EncryptedSubmission.get_encrypted_file_path(i)
        assert enc_path == Path(input_path + ".c4gh")
        enc_path = EncryptedSubmission.get_encryption_header_path(input_path)
        assert enc_path == Path(input_path + ".c4gh_header")
