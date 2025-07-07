import importlib.resources
import itertools
import json
from contextlib import nullcontext
from datetime import date

import pytest
from grz_pydantic_models.mii.consent import Consent
from grz_pydantic_models.submission.metadata.v1 import (
    File,
    FileType,
    GrzSubmissionMetadata,
    ResearchConsent,
    Submission,
)
from pydantic import ValidationError

from . import resources


@pytest.mark.parametrize(
    "dataset,version", itertools.product(["panel", "wes_tumor_germline", "wgs_tumor_germline"], ["1.1.1", "1.1.4"])
)
def test_examples(dataset: str, version: str):
    metadata_str = (
        importlib.resources.files(resources).joinpath("example_metadata", dataset, f"v{version}.json").read_text()
    )
    GrzSubmissionMetadata.model_validate_json(metadata_str)


def test_wgs_trio():
    metadata_str = (
        importlib.resources.files(resources).joinpath("example_metadata", "wgs_trio", "v1.1.7.json").read_text()
    )
    GrzSubmissionMetadata.model_validate_json(metadata_str)


def test_wgs_trio_special_consent():
    """
    Broad Consent obtained before 2025-06-15 for non-index donors is allowed to stand in for mvConsent if missing
    """
    metadata_str = (
        importlib.resources.files(resources)
        .joinpath("example_metadata", "wgs_trio", "v1.1.7.earlyBCException.json")
        .read_text()
    )
    GrzSubmissionMetadata.model_validate_json(metadata_str)

    # only non-index donors can have the special researchConsent exemption
    metadata = json.loads(metadata_str)
    metadata["donors"][0]["mvConsent"]["scope"] = []
    metadata["donors"][0]["researchConsents"][0]["scope"] = metadata["donors"][1]["researchConsents"][0]["scope"]

    with pytest.raises(ValidationError, match="Index donors must have at least a permit of mvSequencing"):
        GrzSubmissionMetadata.model_validate_json(json.dumps(metadata))


def test_example_wgs_lr():
    metadata_str = (
        importlib.resources.files(resources).joinpath("example_metadata", "wgs_lr", "v1.1.4.json").read_text()
    )
    GrzSubmissionMetadata.model_validate_json(metadata_str)


def test_invalid_short_read_submission_with_bam():
    """BAM files should only be allowed in *_lr lab data"""
    metadata = json.loads(
        importlib.resources.files(resources)
        .joinpath("example_metadata", "wgs_tumor_germline", "v1.1.4.json")
        .read_text()
    )
    # add a BAM file
    metadata["donors"][0]["labData"][0]["sequenceData"]["files"].append(
        {
            "filePath": "donor_001/HV5TMDSX7-1-IDUDI0034_S1_L001_R1_001.bam",
            "fileType": "bam",
            "checksumType": "sha256",
            "fileChecksum": "9e87eabc18b726a94a3ffbd8d84df662388bec07b8e3d501ee6a43309c6d43fd",
            "fileSizeInBytes": 129174728987,
            "readLength": 151,
        }
    )

    with pytest.raises(ValidationError):
        GrzSubmissionMetadata.model_validate_json(json.dumps(metadata))


def test_file_extensions():
    File(
        filePath="test/valid.bam",
        fileType=FileType.bam,
        fileChecksum="29647ae83ccac69f2bf4e0f8f37d8f86ad56c578c14432b7a497481031db25b8",
        fileSizeInBytes=0,
        readLength=100,
    )

    with pytest.raises(ValidationError):
        File(
            filePath="test/invalid.bam.gz",
            fileType=FileType.bam,
            fileChecksum="29647ae83ccac69f2bf4e0f8f37d8f86ad56c578c14432b7a497481031db25b8",
            fileSizeInBytes=0,
            readLength=100,
        )


@pytest.mark.parametrize(
    "case,valid",
    (
        ("minimal_consented", True),
        ("minimal_nonconsented", True),
        ("minimal_consented_expired", True),
        ("mii_ig_consent_v2025_example1", True),
        ("invalid_missing_fields", False),
    ),
)
def test_research_consent_parse(case: str, valid: bool):
    expectation = nullcontext() if valid else pytest.raises(ValidationError)

    with expectation:
        Consent.model_validate_json(
            importlib.resources.files(resources).joinpath("example_research_consent", f"{case}.json").read_text()
        )


@pytest.mark.parametrize(
    "cases,consenting",
    (
        (["minimal_consented"], True),
        (["minimal_nonconsented"], False),
        (["minimal_consented_expired"], False),
        (["mii_ig_consent_v2025_example1"], False),
        (["minimal_consented", "minimal_nonconsented"], False),
        (["minimal_consented", "minimal_consented_expired"], True),
        (["minimal_consented", "mii_ig_consent_v2025_example1"], True),
        (["minimal_consented_expired", "mii_ig_consent_v2025_example1"], False),
    ),
)
def test_multi_research_consent(cases: list[str], consenting: bool):
    consents = []
    for case in cases:
        consent = Consent.model_validate_json(
            importlib.resources.files(resources).joinpath("example_research_consent", f"{case}.json").read_text()
        )
        consents.append(ResearchConsent(schemaVersion="2025.0.1", scope=consent))

    assert ResearchConsent.consents_to_research(consents, date=date(year=2025, month=6, day=25)) == consenting


def test_test_submission():
    valid_test_submission = Submission(
        submissionDate="2001-01-01",
        submissionType="initial",
        tanG="TESTGRZXYZ123123456789a6a7cd7ae89d4a9baebdb2a13f85ae9f2a5a62a025",
        localCaseId="CASE12345",
        coverageType="UNK",
        submitterId="123456789",
        genomicDataCenterId="GRZXYZ123",
        clinicalDataNodeId="KDKXYZ123",
        diseaseType="oncological",
        genomicStudyType="single",
        genomicStudySubtype="tumor-only",
        labName="Umbrella",
    )

    # try an invalid test submission date
    submission_date_orig = valid_test_submission.submission_date
    with pytest.raises(ValidationError, match="Test submissions must have a submission date of 01.01.2001"):
        valid_test_submission.submission_date = date(year=2025, month=7, day=1)
    valid_test_submission.submission_date = submission_date_orig

    # try a non-matching GRZ ID
    grz_id_orig = valid_test_submission.genomic_data_center_id
    with pytest.raises(
        ValidationError, match="GRZ ID embedded within test TAN must match the GRZ ID in the submission metadata"
    ):
        valid_test_submission.genomic_data_center_id = "GRZXYZ124"
    valid_test_submission.genomic_data_center_id = grz_id_orig

    # try a non-matching submitter ID
    submitter_id_orig = valid_test_submission.submitter_id
    with pytest.raises(
        ValidationError,
        match="Submitter ID embedded within test TAN must match submitter ID in the submission metadata",
    ):
        valid_test_submission.submitter_id = "123456780"
    valid_test_submission.submitter_id = submitter_id_orig
