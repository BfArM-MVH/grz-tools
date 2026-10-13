"""Tests for the parser module."""

import json
from pathlib import Path

import grz_common.exceptions as grzexc
import pytest
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


@pytest.mark.parametrize(
    "content",
    [b"[]", b'"x"', '{"note": "Grüße"}'.encode("latin-1")],
    ids=["array", "string", "latin-1"],
)
def test_submission_metadata_reports_a_file_that_holds_no_json_object_as_a_validation_error(tmp_path, content):
    """A metadata.json that is valid JSON but no object, or that is not UTF-8, is the submitter's fault."""
    metadata_file = tmp_path / "metadata.json"
    metadata_file.write_bytes(content)

    with pytest.raises(grzexc.SubmissionValidationError):
        SubmissionMetadata(metadata_file)


def test_encrypted_submission():
    input_path = "/submission/files/a.fastq"

    for i in (input_path, Path(input_path)):
        enc_path = EncryptedSubmission.get_encrypted_file_path(i)
        assert enc_path == Path(input_path + ".c4gh")
        enc_path = EncryptedSubmission.get_encryption_header_path(input_path)
        assert enc_path == Path(input_path + ".c4gh_header")


@pytest.mark.parametrize(
    ("metadata_path", "expected_valid"),
    [
        ("tests/mock_files/submissions/valid_submission/metadata/metadata.json", True),
        ("tests/mock_files/no_localCaseId_metadata.json", True),
        ("tests/mock_files/invalid_localCaseId_in_other_field_metadata.json", False),
        ("tests/mock_files/invalid_localCaseId_in_filePath_metadata.json", False),
    ],
    ids=["valid", "no_local_case_id", "local_case_id_in_other", "local_case_id_in_filePath"],
)
def test_local_case_id_in_metadata(metadata_path, expected_valid, identifiers_config_model):
    # Collect errors
    errors = list(SubmissionMetadata(metadata_path).validate(identifiers_config_model.identifiers))

    # Parse errors
    if expected_valid:
        assert [e for e in errors if "localCaseId" in e] == [], errors
    elif metadata_path == "tests/mock_files/invalid_localCaseId_in_other_field_metadata.json":
        assert [
            e
            for e in errors
            if "is the same as localCaseId, which must not appear outside of submission.localCaseId" in e
        ], errors
    else:
        assert [
            e
            for e in errors
            if "contains the localCaseId, which must not appear outside of submission.localCaseId" in e
        ], errors
