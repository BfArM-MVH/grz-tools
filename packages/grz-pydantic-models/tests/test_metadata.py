import copy
import importlib.resources
import itertools
import json
import re
from contextlib import nullcontext
from datetime import date

import pytest
from grz_pydantic_models.mii.consent import Consent
from grz_pydantic_models.submission.metadata import (
    DiseaseType,
    ResearchConsentNoScopeJustification,
)
from grz_pydantic_models.submission.metadata.v1 import (
    File,
    FileType,
    GrzSubmissionMetadata,
    LibraryType,
    ResearchConsent,
)
from grz_pydantic_models_testing import example_metadata, example_research_consent
from packaging.version import Version
from pydantic import ValidationError

TESTED_VERSIONS = ["1.2.1", "1.3.0"]


@pytest.mark.parametrize(
    "dataset,version",
    itertools.product(
        [
            "oncomine_panel_tumor_only",
            "panel_tumor_only",
            "wes_tumor_germline",
            "wgs_tumor_germline",
            "wgs_lr_tumor_only",
            "wgs_trio",
        ],
        TESTED_VERSIONS,
    ),
)
def test_examples(dataset: str, version: str):
    metadata_str = importlib.resources.files(example_metadata).joinpath(dataset, f"v{version}.json").read_text()
    GrzSubmissionMetadata.model_validate_json(metadata_str)


def test_wgs_trio_special_consent():
    """
    Broad Consent obtained before 2025-06-15 for non-index donors is allowed to stand in for mvConsent if missing
    """
    metadata_str = (
        importlib.resources.files(example_metadata).joinpath("wgs_trio", "v1.1.7.earlyBCException.json").read_text()
    )
    GrzSubmissionMetadata.model_validate_json(metadata_str)

    # only non-index donors can have the special researchConsent exemption
    metadata = json.loads(metadata_str)
    metadata["donors"][0]["mvConsent"]["scope"] = []
    metadata["donors"][0]["researchConsents"][0]["scope"] = metadata["donors"][1]["researchConsents"][0]["scope"]

    with pytest.raises(
        ValidationError, match=r"All donors must consent to model project participation for initial submissions."
    ):
        GrzSubmissionMetadata.model_validate_json(json.dumps(metadata))


@pytest.mark.parametrize(
    "version",
    [v for v in TESTED_VERSIONS if Version(v) >= Version("1.1.7")],
)
def test_wgs_trio_no_vcf(version):
    """
    VCFs were downgraded from required to recommended for all submissions.
    """
    metadata_str = importlib.resources.files(example_metadata).joinpath("wgs_trio", f"v{version}.json").read_text()
    GrzSubmissionMetadata.model_validate_json(metadata_str)

    metadata = json.loads(metadata_str)
    # delete the VCF file for the index donor
    del metadata["donors"][0]["labData"][0]["sequenceData"]["files"][2]

    GrzSubmissionMetadata.model_validate_json(json.dumps(metadata))


@pytest.mark.parametrize(
    "version",
    TESTED_VERSIONS,
)
def test_wgs_tumor_germline_missing_dna(version):
    """
    Ensure that both tumor and germline DNA lab data are required for WGS tumor-germline submissions.
    """
    metadata_str = (
        importlib.resources.files(example_metadata).joinpath("wgs_tumor_germline", f"v{version}.json").read_text()
    )

    ### tumor+germline

    # delete germline DNA
    metadata = json.loads(metadata_str)
    assert metadata["donors"][0]["labData"][0]["sequenceSubtype"] == "germline"
    del metadata["donors"][0]["labData"][0]

    # missing germline DNA should fail
    with pytest.raises(
        ValidationError,
        match=r"""Index donor is missing sequence subtypes for submission type 'tumor\+germline': germline""",
    ):
        GrzSubmissionMetadata.model_validate_json(json.dumps(metadata))

    # delete tumor DNA
    metadata = json.loads(metadata_str)
    assert metadata["donors"][0]["labData"][1]["sequenceSubtype"] == "somatic"
    del metadata["donors"][0]["labData"][1]

    # missing tumor DNA should fail
    with pytest.raises(
        ValidationError,
        match=r"""Index donor is missing sequence subtypes for submission type 'tumor\+germline': somatic""",
    ):
        GrzSubmissionMetadata.model_validate_json(json.dumps(metadata))

    ### tumor-only

    # delete germline DNA
    metadata = json.loads(metadata_str)
    metadata["submission"]["genomicStudySubtype"] = "tumor-only"
    assert metadata["donors"][0]["labData"][0]["sequenceSubtype"] == "germline"
    del metadata["donors"][0]["labData"][0]

    # missing germline DNA should pass
    GrzSubmissionMetadata.model_validate_json(json.dumps(metadata))

    # delete tumor DNA
    metadata = json.loads(metadata_str)
    metadata["submission"]["genomicStudySubtype"] = "tumor-only"
    assert metadata["donors"][0]["labData"][1]["sequenceSubtype"] == "somatic"
    del metadata["donors"][0]["labData"][1]

    # missing tumor DNA should fail
    with pytest.raises(
        ValidationError,
        match=r"""Index donor is missing sequence subtypes for submission type 'tumor-only': somatic""",
    ):
        GrzSubmissionMetadata.model_validate_json(json.dumps(metadata))

    ### germline-only

    # delete germline DNA
    metadata = json.loads(metadata_str)
    metadata["submission"]["genomicStudySubtype"] = "germline-only"
    assert metadata["donors"][0]["labData"][0]["sequenceSubtype"] == "germline"
    del metadata["donors"][0]["labData"][0]

    # missing germline DNA should fail
    with pytest.raises(
        ValidationError,
        match=r"""Index donor is missing sequence subtypes for submission type 'germline-only': germline""",
    ):
        GrzSubmissionMetadata.model_validate_json(json.dumps(metadata))

    # delete tumor DNA
    metadata = json.loads(metadata_str)
    metadata["submission"]["genomicStudySubtype"] = "germline-only"
    assert metadata["donors"][0]["labData"][1]["sequenceSubtype"] == "somatic"
    del metadata["donors"][0]["labData"][1]

    # missing tumor DNA should pass
    GrzSubmissionMetadata.model_validate_json(json.dumps(metadata))


@pytest.mark.parametrize(
    "version",
    [v for v in TESTED_VERSIONS if Version(v) >= Version("1.3.0")],
)
def test_wgs_trio_1_3_fail_empty_consent_list(version: str):
    """As of v1.3, empty consent lists are no longer allowed."""
    metadata_str = importlib.resources.files(example_metadata).joinpath("wgs_trio", f"v{version}.json").read_text()
    metadata = json.loads(metadata_str)
    metadata["donors"][0]["researchConsents"] = []
    with pytest.raises(ValidationError):
        GrzSubmissionMetadata.model_validate_json(metadata)


@pytest.mark.parametrize(
    "version",
    [v for v in TESTED_VERSIONS if Version(v) >= Version("1.3.0")],
)
def test_wgs_trio_1_3_fail_malformed_consent(version: str):
    """As of v1.3, non-empty scope or noScopeJustification must be provided."""
    metadata_str = importlib.resources.files(example_metadata).joinpath("wgs_trio", f"v{version}.json").read_text()
    metadata = json.loads(metadata_str)

    # scope, if provided, must be a valid consent object
    del metadata["donors"][0]["researchConsents"][0]["scope"]["scope"]
    with pytest.raises(ValidationError, match=r"scope must be a valid MII Broad Consent as of metadata v1.3"):
        GrzSubmissionMetadata.model_validate_json(json.dumps(metadata))

    # scope can't be an empty dict
    metadata["donors"][0]["researchConsents"][0]["scope"] = {}
    with pytest.raises(ValidationError):
        GrzSubmissionMetadata.model_validate_json(json.dumps(metadata))

    metadata["donors"][0]["researchConsents"][0]["noScopeJustification"] = "patient unable to consent"
    # scope, even if an empty dict, can't be provided along with noScopeJustification
    with pytest.raises(ValidationError):
        GrzSubmissionMetadata.model_validate_json(json.dumps(metadata))

    del metadata["donors"][0]["researchConsents"][0]["scope"]
    GrzSubmissionMetadata.model_validate_json(json.dumps(metadata))

    # schemaVersion can be missing now
    del metadata["donors"][0]["researchConsents"][0]["schemaVersion"]
    GrzSubmissionMetadata.model_validate_json(json.dumps(metadata))

    # but presentationDate no longer can
    del metadata["donors"][0]["researchConsents"][0]["presentationDate"]
    with pytest.raises(ValidationError):
        GrzSubmissionMetadata.model_validate_json(json.dumps(metadata))


@pytest.mark.parametrize(
    "dataset,version",
    itertools.product(["panel_tumor_only", "wes_tumor_germline", "wgs_tumor_germline", "wgs_trio"], TESTED_VERSIONS),
)
def test_invalid_short_read_submission_with_bam(dataset: str, version: str):
    """BAM files should only be allowed in *_lr lab data"""
    metadata = json.loads(importlib.resources.files(example_metadata).joinpath(dataset, f"v{version}.json").read_text())
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


@pytest.mark.parametrize("version", TESTED_VERSIONS)
def test_index_rna_without_dna(version: str):
    """Donors can only have RNA data if DNA data also present."""
    metadata = json.loads(
        importlib.resources.files(example_metadata).joinpath("wes_tumor_germline", f"v{version}.json").read_text()
    )
    # reduce to a single lab datum
    metadata["donors"][0]["labData"] = [metadata["donors"][0]["labData"][0]]
    # set the library type to RNA
    metadata["donors"][0]["labData"][0]["libraryType"] = "wxs"
    metadata["donors"][0]["labData"][0]["sequenceType"] = "rna"

    with pytest.raises(
        ValidationError, match="Index donor must have at least one lab datum with one of the following library types"
    ):
        GrzSubmissionMetadata.model_validate_json(json.dumps(metadata))


@pytest.mark.parametrize("version", TESTED_VERSIONS)
def test_index_rna_with_dna(version: str):
    """Donors can only have RNA data if DNA data also present."""
    metadata = json.loads(
        importlib.resources.files(example_metadata).joinpath("wes_tumor_germline", f"v{version}.json").read_text()
    )
    # duplicate the last lab datum
    metadata["donors"][0]["labData"].append(copy.deepcopy(metadata["donors"][0]["labData"][-1]))
    # set the library type to RNA
    metadata["donors"][0]["labData"][-1]["libraryType"] = "wxs"
    metadata["donors"][0]["labData"][-1]["sequenceType"] = "rna"
    metadata["donors"][0]["labData"][-1]["labDataName"] = metadata["donors"][0]["labData"][-2]["labDataName"] + " RNA"

    # fix file checksums to be different
    import hashlib

    for file_idx, file in enumerate(metadata["donors"][0]["labData"][-1]["sequenceData"]["files"]):
        file["fileChecksum"] = hashlib.sha256(hex(file_idx).encode("utf8")).hexdigest()

    GrzSubmissionMetadata.model_validate_json(json.dumps(metadata))


@pytest.mark.parametrize("version", TESTED_VERSIONS)
def test_lab_datum(version: str):
    metadata = GrzSubmissionMetadata.model_validate_json(
        importlib.resources.files(example_metadata).joinpath("wes_tumor_germline", f"v{version}.json").read_text()
    )
    with pytest.raises(ValueError, match=re.escape("Long read libraries can't be paired-end.")):
        metadata.donors[0].lab_data[0].library_type = "wes_lr"


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
        ("minimal_consented_with_datetime", True),
        ("minimal_consented_with_nonzero_datetime", True),
        ("extra_consented", True),
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
            importlib.resources.files(example_research_consent).joinpath(f"{case}.json").read_text()
        )


def test_research_consent_tolerates_fhir_extensions():
    """FHIR allows id/extension on every element, so they must not be rejected on codings, concepts, or periods."""
    consent_raw = json.loads(
        importlib.resources.files(example_research_consent).joinpath("mii_ig_consent_v2025_example1.json").read_text()
    )

    consent_raw["scope"]["coding"][0]["extension"] = [{"url": "http://example.org/ext", "valueString": "x"}]
    consent_raw["scope"]["coding"][0]["id"] = "coding-1"
    consent_raw["category"][0]["extension"] = [{"url": "http://example.org/ext"}]
    consent_raw["provision"]["provision"][0]["period"]["extension"] = [{"url": "http://example.org/ext"}]
    consent_raw["provision"]["provision"][0]["period"]["id"] = "period-1"

    Consent.model_validate(consent_raw)


def test_research_consent_still_rejects_unknown_fields():
    """Allowing id/extension must not turn into accepting arbitrary unknown fields."""
    consent_raw = json.loads(
        importlib.resources.files(example_research_consent).joinpath("mii_ig_consent_v2025_example1.json").read_text()
    )
    consent_raw["scope"]["coding"][0]["pandorras-box"] = "I am in"

    with pytest.raises(ValidationError):
        Consent.model_validate(consent_raw)


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
            importlib.resources.files(example_research_consent).joinpath(f"{case}.json").read_text()
        )
        consents.append(ResearchConsent(schemaVersion="2025.0.1", scope=consent))

    assert ResearchConsent.consents_to_research(consents, date=date(year=2025, month=6, day=25)) == consenting


def test_research_consent_subprovisions_deny_permit():
    """Within one research consent's subprovisions, deny before permit should return a non-consented state."""
    consent_raw = json.loads(
        importlib.resources.files(example_research_consent).joinpath("minimal_nonconsented.json").read_text()
    )

    # add a permit subprovision object for same consent object, after the deny subprovision
    new_permit_subprovision = copy.deepcopy(consent_raw["provision"]["provision"][0])
    new_permit_subprovision["type"] = "permit"
    consent_raw["provision"]["provision"].append(new_permit_subprovision)

    consent = Consent.model_validate_json(json.dumps(consent_raw))

    assert not ResearchConsent.consents_to_research(
        [ResearchConsent(scope=consent)], date=date(year=2025, month=10, day=13)
    )


def test_research_consents_deny_permit():
    """Having two research consents, where deny comes before permit, should return a non-consented state."""
    consent_raw = json.loads(
        importlib.resources.files(example_research_consent).joinpath("minimal_nonconsented.json").read_text()
    )
    consent1 = Consent.model_validate_json(json.dumps(consent_raw))

    # add a permit consent object for same donor
    consent_raw["provision"]["provision"][0]["type"] = "permit"
    consent2 = Consent.model_validate_json(json.dumps(consent_raw))

    assert not ResearchConsent.consents_to_research(
        (ResearchConsent(scope=consent1), ResearchConsent(scope=consent2)), date=date(year=2025, month=10, day=13)
    )


def test_research_consent_no_subprovisions():
    """Consent objects are allowed to have no provisions under the root."""
    consent_json_raw = json.loads(
        importlib.resources.files(example_research_consent).joinpath("minimal_consented.json").read_text()
    )
    del consent_json_raw["provision"]["provision"]
    Consent.model_validate_json(json.dumps(consent_json_raw))


@pytest.mark.parametrize(
    "dataset,version",
    itertools.product(
        [
            "oncomine_panel_tumor_only",
            "panel_tumor_only",
            "wes_tumor_germline",
        ],
        TESTED_VERSIONS,
    ),
)
def test_disease_type_rare_missing_wgs_raises(dataset: str, version: str):
    """
    These datasets natively use panel or wes. For rare diseases, they fail the wgs check.
    """
    metadata = json.loads(importlib.resources.files(example_metadata).joinpath(dataset, f"v{version}.json").read_text())
    metadata["submission"]["diseaseType"] = DiseaseType.rare.value
    metadata["submission"]["submissionDate"] = "2026-06-01"

    with pytest.raises(
        ValidationError, match=r"the index donor must have at least one lab datum with libraryType 'wgs' or 'wgs_lr'"
    ):
        GrzSubmissionMetadata.model_validate_json(json.dumps(metadata))


@pytest.mark.parametrize(
    "dataset,version",
    itertools.product(["wes_tumor_germline"], TESTED_VERSIONS),
)
def test_disease_type_rare_missing_wgs_warns_before_cutoff(dataset: str, version: str, caplog):
    """Before the cutoff, missing wgs for a rare disease should only log a warning."""
    metadata = json.loads(importlib.resources.files(example_metadata).joinpath(dataset, f"v{version}.json").read_text())
    metadata["submission"]["diseaseType"] = DiseaseType.rare.value
    metadata["submission"]["submissionDate"] = "2026-05-31"

    try:
        GrzSubmissionMetadata.model_validate_json(json.dumps(metadata))
    except ValidationError as e:
        assert "starting 01.06.2026" not in str(e)

    assert "starting 01.06.2026" in caplog.text


@pytest.mark.parametrize("version", TESTED_VERSIONS)
def test_disease_type_rare_index_has_wgs_and_panel_passes(version: str):
    """
    If the index donor has at least WGS but additionally some panel data, it should still pass.
    """
    metadata = json.loads(
        importlib.resources.files(example_metadata).joinpath("wgs_trio", f"v{version}.json").read_text()
    )
    metadata["submission"]["diseaseType"] = DiseaseType.rare.value
    metadata["submission"]["submissionDate"] = "2026-06-01"

    # duplicate the WGS lab datum but change its type to panel
    panel_datum = copy.deepcopy(metadata["donors"][0]["labData"][0])
    panel_datum["libraryType"] = LibraryType.panel.value
    panel_datum["labDataName"] += "_panel"

    # panels require a BED file
    panel_datum["sequenceData"]["files"].append(
        {"filePath": "dummy_panel_target.bed", "fileType": "bed", "fileChecksum": "c" * 64, "fileSizeInBytes": 1000}
    )

    # change file paths/checksums so it doesn't fail unique run/checksum validators
    for i, file in enumerate(panel_datum["sequenceData"]["files"]):
        file["fileChecksum"] = f"c{i:063d}"
        if file["fileType"] == "fastq":
            file["filePath"] = f"panel_{i}.fastq.gz"

    metadata["donors"][0]["labData"].append(panel_datum)

    try:
        GrzSubmissionMetadata.model_validate_json(json.dumps(metadata))
    except ValidationError as e:
        assert "starting 01.06.2026" not in str(e)


@pytest.mark.parametrize(
    "dataset,version",
    itertools.product(
        [
            "wgs_tumor_germline",
            "wgs_lr_tumor_only",
            "wgs_trio",
        ],
        TESTED_VERSIONS,
    ),
)
def test_disease_type_rare_valid_library_passes(dataset: str, version: str):
    """These datasets natively use wgs or wgs_lr. For rare diseases, they must pass."""
    metadata = json.loads(importlib.resources.files(example_metadata).joinpath(dataset, f"v{version}.json").read_text())
    metadata["submission"]["diseaseType"] = DiseaseType.rare.value
    metadata["submission"]["submissionDate"] = "2026-06-01"

    try:
        GrzSubmissionMetadata.model_validate_json(json.dumps(metadata))
    except ValidationError as e:
        assert "starting 01.06.2026" not in str(e)


@pytest.mark.parametrize("version", TESTED_VERSIONS)
def test_disease_type_rare_non_index_unrestricted(version: str):
    """Non-index donors (e.g., parents in a trio) are exempt from the WGS rule."""
    metadata = json.loads(
        importlib.resources.files(example_metadata).joinpath("wgs_trio", f"v{version}.json").read_text()
    )
    metadata["submission"]["diseaseType"] = DiseaseType.rare.value
    metadata["submission"]["submissionDate"] = "2026-06-01"

    # modify a non-index donor (donor[1]) to use a panel
    for idx, datum in enumerate(metadata["donors"][1]["labData"]):
        datum["libraryType"] = LibraryType.panel.value

        # panels require a BED file
        datum["sequenceData"]["files"].append(
            {
                "filePath": f"dummy_target_{idx}.bed",
                "fileType": "bed",
                "fileChecksum": f"c{idx:063d}",
                "fileSizeInBytes": 1000,
            }
        )

    try:
        # should pass because the index donor still has WGS
        GrzSubmissionMetadata.model_validate_json(json.dumps(metadata))
    except ValidationError as e:
        assert "starting 01.06.2026" not in str(e)


@pytest.mark.parametrize(
    "version,justification",
    itertools.product(
        TESTED_VERSIONS,
        [ResearchConsentNoScopeJustification.LE_TECH.value, ResearchConsentNoScopeJustification.LE_ORG.value],
    ),
)
def test_no_scope_justification_tech_org_raises(version: str, justification: str):
    """LE_TECH and LE_ORG justifications should raise an error starting exactly 01.06.2026."""
    metadata = json.loads(
        importlib.resources.files(example_metadata).joinpath("wgs_trio", f"v{version}.json").read_text()
    )
    metadata["submission"]["submissionDate"] = "2026-06-01"

    # Apply to all consents to be thorough
    for donor in metadata["donors"]:
        for consent in donor.get("researchConsents", []):
            consent["noScopeJustification"] = justification
            consent.pop("scope", None)

    with pytest.raises(ValidationError, match=r"is no longer allowed starting 01\.06\.2026"):
        GrzSubmissionMetadata.model_validate_json(json.dumps(metadata))


@pytest.mark.parametrize(
    "version,justification",
    itertools.product(
        TESTED_VERSIONS,
        [ResearchConsentNoScopeJustification.LE_TECH.value, ResearchConsentNoScopeJustification.LE_ORG.value],
    ),
)
def test_no_scope_justification_tech_org_warns(version: str, justification: str, caplog):
    """LE_TECH and LE_ORG justifications should only warn strictly before 01.06.2026."""
    metadata = json.loads(
        importlib.resources.files(example_metadata).joinpath("wgs_trio", f"v{version}.json").read_text()
    )
    metadata["submission"]["submissionDate"] = "2026-05-31"

    for donor in metadata["donors"]:
        for consent in donor.get("researchConsents", []):
            consent["noScopeJustification"] = justification
            consent.pop("scope", None)

    try:
        GrzSubmissionMetadata.model_validate_json(json.dumps(metadata))
    except ValidationError as e:
        assert "is no longer allowed starting 01.06.2026" not in str(e)

    assert "is no longer allowed starting 01.06.2026" in caplog.text


@pytest.mark.parametrize(
    "version,justification",
    itertools.product(
        TESTED_VERSIONS,
        [
            ResearchConsentNoScopeJustification.UNABLE.value,
            ResearchConsentNoScopeJustification.REFUSED.value,
            ResearchConsentNoScopeJustification.NO_RETURN.value,
            ResearchConsentNoScopeJustification.OTHER.value,
        ],
    ),
)
def test_no_scope_justification_standard_passes_after_cutoff(version: str, justification: str):
    """Standard valid justifications should continue to pass seamlessly after 01.06.2026."""
    metadata = json.loads(
        importlib.resources.files(example_metadata).joinpath("wgs_trio", f"v{version}.json").read_text()
    )
    metadata["submission"]["submissionDate"] = "2026-06-01"

    for donor in metadata["donors"]:
        for consent in donor.get("researchConsents", []):
            consent["noScopeJustification"] = justification
            consent.pop("scope", None)

    try:
        GrzSubmissionMetadata.model_validate_json(json.dumps(metadata))
    except ValidationError as e:
        assert "is no longer allowed starting 01.06.2026" not in str(e)
