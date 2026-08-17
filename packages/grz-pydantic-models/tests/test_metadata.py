import copy
import importlib.resources
import itertools
import json
import re
from datetime import UTC, date, datetime, timedelta, timezone

import pytest
from grz_pydantic_models.mii.consent import (
    BROAD_CONSENT_DOCUMENT_OIDS,
    EXPECTED_CATEGORIES,
    EXPECTED_SCOPE_CODING_CODE,
    EXPECTED_SCOPE_CODING_SYSTEM,
    MII_BROAD_CONSENT_OID,
    MII_CONSENT_CATEGORY_SYSTEMS,
    MII_CONSENT_VERSION_MODULES_SYSTEM,
    BroadConsentVersion,
    Consent,
    ConsentDocumentKind,
    ConsentProvision,
    FhirDateTime,
    Identifier,
    Period,
    RootConsentProvision,
)
from grz_pydantic_models.submission.metadata import (
    DiseaseType,
    ResearchConsentNoScopeJustification,
    get_accepted_versions,
)
from grz_pydantic_models.submission.metadata.v1 import (
    PROFILES_REQUIRING_PERIOD_END,
    RESEARCH_CONSENT_PACKAGE_PROFILES,
    RESEARCH_CONSENT_SCHEMA_VERSIONS,
    File,
    FileType,
    GrzSubmissionMetadata,
    LibraryType,
    ResearchConsent,
    ResearchConsentCodes,
)
from grz_pydantic_models_testing import example_metadata, example_research_consent, example_terminology
from packaging.version import Version
from pydantic import ValidationError

TESTED_VERSIONS = ["1.2.1", "1.3.0"]


def _metadata_raw(dataset: str, version: str) -> str:
    """The raw JSON text of an example metadata submission."""
    return importlib.resources.files(example_metadata).joinpath(dataset, f"v{version}.json").read_text()


def _metadata(dataset: str, version: str) -> GrzSubmissionMetadata:
    """An example metadata submission parsed as the model under test."""
    return GrzSubmissionMetadata.model_validate_json(_metadata_raw(dataset, version))


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
    _metadata(dataset, version)


def test_wgs_trio_special_consent():
    """
    Broad Consent obtained before 2025-06-15 for non-index donors is allowed to stand in for mvConsent if missing
    """
    metadata_str = _metadata_raw("wgs_trio", "1.1.7.earlyBCException")
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
    metadata_str = _metadata_raw("wgs_trio", version)
    GrzSubmissionMetadata.model_validate_json(metadata_str)

    metadata = json.loads(metadata_str)
    # delete the VCF file for the index donor
    del metadata["donors"][0]["labData"][0]["sequenceData"]["files"][2]

    GrzSubmissionMetadata.model_validate_json(json.dumps(metadata))


@pytest.mark.parametrize(
    "genomic_study_subtype,deleted_subtype,should_raise",
    (
        ("tumor+germline", "germline", True),
        ("tumor+germline", "somatic", True),
        ("tumor-only", "germline", False),
        ("tumor-only", "somatic", True),
        ("germline-only", "germline", True),
        ("germline-only", "somatic", False),
    ),
)
@pytest.mark.parametrize("version", TESTED_VERSIONS)
def test_wgs_tumor_germline_missing_dna(
    version: str, genomic_study_subtype: str, deleted_subtype: str, should_raise: bool
):
    """
    Ensure that both tumor and germline DNA lab data are required for WGS tumor-germline submissions,
    except where the deleted subtype is not required by the declared genomicStudySubtype.
    """
    metadata = json.loads(_metadata_raw("wgs_tumor_germline", version))
    # The example predates the 04.12.2025 cutoff, before which a missing subtype only warns.
    metadata["submission"]["submissionDate"] = "2025-12-04"
    metadata["submission"]["genomicStudySubtype"] = genomic_study_subtype

    lab_datum_index = {"germline": 0, "somatic": 1}[deleted_subtype]
    assert metadata["donors"][0]["labData"][lab_datum_index]["sequenceSubtype"] == deleted_subtype
    del metadata["donors"][0]["labData"][lab_datum_index]

    if should_raise:
        with pytest.raises(
            ValidationError,
            match=rf"""Index donor is missing sequence subtypes for submission type '{re.escape(genomic_study_subtype)}': {deleted_subtype}""",
        ):
            GrzSubmissionMetadata.model_validate_json(json.dumps(metadata))
    else:
        GrzSubmissionMetadata.model_validate_json(json.dumps(metadata))


@pytest.mark.parametrize("version", TESTED_VERSIONS)
def test_missing_sequence_subtype_warns_before_cutoff(version: str, caplog):
    """
    The labData/genomicStudySubtype check was published in grz-pydantic-models v2.5.0 (04.12.2025).
    Submissions predating that release must only warn, so that backfills and other operations
    replaying historical data do not fail retroactively.
    """
    metadata = json.loads(_metadata_raw("wgs_tumor_germline", version))
    metadata["submission"]["submissionDate"] = "2025-12-03"
    assert metadata["donors"][0]["labData"][0]["sequenceSubtype"] == "germline"
    del metadata["donors"][0]["labData"][0]

    # must not raise for this rule
    GrzSubmissionMetadata.model_validate_json(json.dumps(metadata))

    assert "Index donor is missing sequence subtypes" in caplog.text
    assert "starting 04.12.2025" in caplog.text


@pytest.mark.parametrize("version", TESTED_VERSIONS)
def test_missing_sequence_subtype_raises_on_cutoff(version: str):
    """On the cutoff date itself the check is enforced."""
    metadata = json.loads(_metadata_raw("wgs_tumor_germline", version))
    metadata["submission"]["submissionDate"] = "2025-12-04"
    assert metadata["donors"][0]["labData"][0]["sequenceSubtype"] == "germline"
    del metadata["donors"][0]["labData"][0]

    with pytest.raises(
        ValidationError,
        match=r"""Index donor is missing sequence subtypes for submission type 'tumor\+germline': germline""",
    ):
        GrzSubmissionMetadata.model_validate_json(json.dumps(metadata))


@pytest.mark.parametrize(
    "version",
    [v for v in TESTED_VERSIONS if Version(v) >= Version("1.3.0")],
)
def test_wgs_trio_1_3_fail_empty_consent_list(version: str):
    """As of v1.3, empty consent lists are no longer allowed."""
    metadata = json.loads(_metadata_raw("wgs_trio", version))
    metadata["donors"][0]["researchConsents"] = []
    with pytest.raises(ValidationError, match=r"Donors must have research consent as of metadata schema v1\.3"):
        GrzSubmissionMetadata.model_validate_json(json.dumps(metadata))


@pytest.mark.parametrize(
    "version",
    [v for v in TESTED_VERSIONS if Version(v) >= Version("1.3.0")],
)
def test_wgs_trio_1_3_fail_scope_must_be_a_consent(version: str):
    """As of v1.3, a scope that fails to parse as a Consent is rejected."""
    metadata = json.loads(_metadata_raw("wgs_trio", version))
    del metadata["donors"][0]["researchConsents"][0]["scope"]["scope"]

    with pytest.raises(ValidationError, match=r"scope must be a valid MII Broad Consent as of metadata v1.3"):
        GrzSubmissionMetadata.model_validate_json(json.dumps(metadata))


@pytest.mark.parametrize(
    "version",
    [v for v in TESTED_VERSIONS if Version(v) >= Version("1.3.0")],
)
def test_wgs_trio_1_3_fail_empty_scope(version: str):
    """As of v1.3, scope can't be an empty dict."""
    metadata = json.loads(_metadata_raw("wgs_trio", version))
    metadata["donors"][0]["researchConsents"][0]["scope"] = {}

    with pytest.raises(ValidationError):
        GrzSubmissionMetadata.model_validate_json(json.dumps(metadata))


@pytest.mark.parametrize(
    "version",
    [v for v in TESTED_VERSIONS if Version(v) >= Version("1.3.0")],
)
def test_wgs_trio_1_3_fail_scope_and_justification_mutually_exclusive(version: str):
    """As of v1.3, scope, even if an empty dict, can't be provided along with noScopeJustification."""
    metadata = json.loads(_metadata_raw("wgs_trio", version))
    metadata["donors"][0]["researchConsents"][0]["scope"] = {}
    metadata["donors"][0]["researchConsents"][0]["noScopeJustification"] = "patient unable to consent"

    with pytest.raises(ValidationError):
        GrzSubmissionMetadata.model_validate_json(json.dumps(metadata))


@pytest.mark.parametrize(
    "version",
    [v for v in TESTED_VERSIONS if Version(v) >= Version("1.3.0")],
)
def test_wgs_trio_1_3_justification_replaces_scope(version: str):
    """As of v1.3, a justification may stand in for a scope, with schemaVersion still present."""
    metadata = json.loads(_metadata_raw("wgs_trio", version))
    metadata["donors"][0]["researchConsents"][0]["noScopeJustification"] = "patient unable to consent"
    del metadata["donors"][0]["researchConsents"][0]["scope"]

    GrzSubmissionMetadata.model_validate_json(json.dumps(metadata))


@pytest.mark.parametrize(
    "version",
    [v for v in TESTED_VERSIONS if Version(v) >= Version("1.3.0")],
)
def test_wgs_trio_1_3_schema_version_now_optional(version: str):
    """As of v1.3, researchConsent no longer requires schemaVersion."""
    metadata = json.loads(_metadata_raw("wgs_trio", version))
    metadata["donors"][0]["researchConsents"][0]["noScopeJustification"] = "patient unable to consent"
    del metadata["donors"][0]["researchConsents"][0]["scope"]
    del metadata["donors"][0]["researchConsents"][0]["schemaVersion"]

    GrzSubmissionMetadata.model_validate_json(json.dumps(metadata))


@pytest.mark.parametrize(
    "version",
    [v for v in TESTED_VERSIONS if Version(v) >= Version("1.3.0")],
)
def test_wgs_trio_1_3_fail_missing_presentation_date(version: str):
    """As of v1.3, presentationDate is required even where schemaVersion is not."""
    metadata = json.loads(_metadata_raw("wgs_trio", version))
    metadata["donors"][0]["researchConsents"][0]["noScopeJustification"] = "patient unable to consent"
    del metadata["donors"][0]["researchConsents"][0]["scope"]
    del metadata["donors"][0]["researchConsents"][0]["schemaVersion"]
    del metadata["donors"][0]["researchConsents"][0]["presentationDate"]

    with pytest.raises(ValidationError):
        GrzSubmissionMetadata.model_validate_json(json.dumps(metadata))


@pytest.mark.parametrize(
    "dataset,version",
    itertools.product(["panel_tumor_only", "wes_tumor_germline", "wgs_tumor_germline", "wgs_trio"], TESTED_VERSIONS),
)
def test_invalid_short_read_submission_with_bam(dataset: str, version: str):
    """BAM files should only be allowed in *_lr lab data"""
    metadata = json.loads(_metadata_raw(dataset, version))
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
    metadata = json.loads(_metadata_raw("wes_tumor_germline", version))
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
    metadata = json.loads(_metadata_raw("wes_tumor_germline", version))
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
    submission = _metadata("wes_tumor_germline", version)
    with pytest.raises(ValueError, match=re.escape("Long read libraries can't be paired-end.")):
        submission.donors[0].lab_data[0].library_type = "wes_lr"


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


def _consent_raw(case: str) -> dict:
    """The raw JSON of an example consent, for tests that mutate it before validation."""
    return json.loads(importlib.resources.files(example_research_consent).joinpath(f"{case}.json").read_text())


def _consent(case: str) -> Consent:
    """An example consent parsed as the model under test."""
    return Consent.model_validate(_consent_raw(case))


# Example consents the MII consent profile permits; each must parse.
VALID_CONSENT_CASES = (
    "minimal_consented",
    "minimal_consented_with_datetime",
    "minimal_consented_with_nonzero_datetime",
    "extra_consented",
    "minimal_nonconsented",
    "minimal_consented_expired",
    "mii_ig_consent_v2025_example1",
    # examples shipped inside the MII consent package 2026.0.0, vendored verbatim
    "mii_ig_consent_v2026_example1",
    "mii_ig_consent_v2026_example2",
    "mii_ig_consent_v2026_example3",
    "mii_ig_consent_v2026_example4",
    # instances the MII consent profile permits and that must therefore not be rejected
    "minimal_consented_open_ended",
    "minimal_consented_version_modules_category",
    "minimal_consented_extra_multi_coding_category",
    "minimal_consented_patient_by_identifier",
    "minimal_consented_bc_v1_7_2",
    "minimal_consented_bc_v1_6d_bare_oid",
    "minimal_consented_unknown_policy_oid",
    "withdrawal_complete",
    "rejection_bc_v1_7_2",
    "status_rejected",
)
# Instances the MII consent profile forbids and whose contents would otherwise be dropped silently.
INVALID_CONSENT_CASES = (
    "invalid_missing_fields",
    "invalid_wrong_resource_type",
    "invalid_third_level_provision",
    "invalid_root_provision_code",
    "invalid_mii_category_multi_coding",
)


def test_consent_case_lists_match_the_fixture_package():
    """Every fixture must be listed and every listed fixture must ship, so drift fails loudly."""
    shipped = {
        entry.name.removesuffix(".json")
        for entry in importlib.resources.files(example_research_consent).iterdir()
        if entry.name.endswith(".json")
    }
    assert shipped == set(VALID_CONSENT_CASES) | set(INVALID_CONSENT_CASES)


@pytest.mark.parametrize("case", VALID_CONSENT_CASES)
def test_research_consent_parses(case: str):
    _consent(case)


@pytest.mark.parametrize("case", INVALID_CONSENT_CASES)
def test_research_consent_rejected(case: str):
    with pytest.raises(ValidationError):
        _consent(case)


def test_research_consent_tolerates_fhir_extensions():
    """FHIR allows id/extension on every element, so they must not be rejected on codings, concepts, or periods."""
    consent_raw = _consent_raw("mii_ig_consent_v2025_example1")

    consent_raw["scope"]["coding"][0]["extension"] = [{"url": "http://example.org/ext", "valueString": "x"}]
    consent_raw["scope"]["coding"][0]["id"] = "coding-1"
    consent_raw["category"][0]["extension"] = [{"url": "http://example.org/ext"}]
    consent_raw["provision"]["provision"][0]["period"]["extension"] = [{"url": "http://example.org/ext"}]
    consent_raw["provision"]["provision"][0]["period"]["id"] = "period-1"

    Consent.model_validate(consent_raw)


def test_research_consent_parses_verification_date():
    """A well-formed FHIR verification block (verified + verificationDate) must be parsed, not ignored."""
    consent_raw = _consent_raw("mii_ig_consent_v2025_example1")
    consent_raw["verification"] = [{"verified": True, "verificationDate": "2026-03-27T11:27:01+01:00"}]

    consent = Consent.model_validate(consent_raw)

    assert consent.verification is not None
    assert consent.verification[0].verified is True
    assert consent.verification[0].verification_date is not None


def test_research_consent_rejects_verification_without_verified():
    """`verified` is required (FHIR 1..1), so a verification entry missing it must be rejected."""
    consent_raw = _consent_raw("mii_ig_consent_v2025_example1")
    consent_raw["verification"] = [{"verificationDate": "2026-03-27T11:27:01+01:00"}]

    with pytest.raises(ValidationError):
        Consent.model_validate(consent_raw)


def test_research_consent_still_rejects_unknown_fields():
    """Allowing id/extension must not turn into accepting arbitrary unknown fields."""
    consent_raw = _consent_raw("mii_ig_consent_v2025_example1")
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
    consents = [ResearchConsent(schemaVersion="2025.0.1", scope=_consent(case)) for case in cases]

    assert ResearchConsent.consents_to_research(consents, date=date(year=2025, month=6, day=25)) == consenting


@pytest.mark.parametrize("version", RESEARCH_CONSENT_SCHEMA_VERSIONS)
def test_every_published_package_version_parses_the_same_consent(version: str):
    """Every version the model claims to support must accept one and the same consent."""
    research_consent = ResearchConsent(schemaVersion=version, scope=_consent("minimal_consented"))
    assert isinstance(research_consent.scope, Consent)
    assert ResearchConsent.consents_to_research([research_consent], date=date(year=2024, month=1, day=1))


@pytest.mark.parametrize(
    "version",
    (
        "1.0.7",  # a consent package predating the ones this model covers
        "9999.0.0",  # the accepted versions are an allow list, not a floor
        "not-a-version",
    ),
)
def test_research_consent_schema_version_rejected(version: str):
    """Anything outside the allow list is refused, including versions newer than any known one."""
    with pytest.raises(ValidationError):
        ResearchConsent(schemaVersion=version, scope=_consent("minimal_consented"))


def test_metadata_declaring_schema_1_3_1_validates():
    """Schema 1.3.1 only widens the consent version enum, so it must validate like 1.3.0."""
    metadata = json.loads(_metadata_raw("wgs_trio", "1.3.0"))
    metadata["$schema"] = metadata["$schema"].replace("/v1.3.0/", "/v1.3.1/")

    submission = GrzSubmissionMetadata.model_validate(metadata)
    assert submission.get_schema_version() == "1.3.1"
    # grz-common's validate() gates on this set, so parsing alone does not make 1.3.1 accepted
    assert submission.get_schema_version() in get_accepted_versions()


def test_research_consent_schema_version_error_names_the_accepted_versions():
    """The message reaches submitters, so it must say what is accepted rather than only what is not."""
    with pytest.raises(ValidationError) as excinfo:
        ResearchConsent(schemaVersion="9999.0.0", scope=_consent("minimal_consented"))

    message = str(excinfo.value)
    for version in RESEARCH_CONSENT_SCHEMA_VERSIONS:
        assert version in message


@pytest.mark.parametrize(
    "version,category_system",
    itertools.product(
        [v for v in TESTED_VERSIONS if Version(v) >= Version("1.3.0")],
        MII_CONSENT_CATEGORY_SYSTEMS,
    ),
)
def test_wgs_trio_accepts_renamed_consent_category_system(version: str, category_system: str):
    """
    Either category CodeSystem spelling must parse: the 2026.0.0 package renamed it but still
    ships examples using the old name.
    """
    metadata = json.loads(_metadata_raw("wgs_trio", version))
    research_consent = metadata["donors"][0]["researchConsents"][0]
    mii_codings = [
        coding
        for category in research_consent["scope"]["category"]
        for coding in category["coding"]
        if coding["code"] == MII_BROAD_CONSENT_OID
    ]
    assert mii_codings, "example consent should contain the MII broad consent category"
    for coding in mii_codings:
        coding["system"] = category_system

    submission = GrzSubmissionMetadata.model_validate_json(json.dumps(metadata))
    assert isinstance(submission.donors[0].research_consents[0].scope, Consent)


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
def test_example_research_consent_scopes_parse_as_consent(dataset: str, version: str):
    """The scope union falls back to dict silently, so every example scope must parse as a Consent."""
    submission = _metadata(dataset, version)
    for donor in submission.donors:
        for research_consent in donor.research_consents:
            assert research_consent.scope is None or isinstance(research_consent.scope, Consent)


def test_research_consent_open_ended_provision_period():
    """A provision without an end must parse and stay in force, not count as expired."""
    consent_raw = _consent_raw("minimal_consented")
    del consent_raw["provision"]["period"]["end"]
    for provision in consent_raw["provision"]["provision"]:
        del provision["period"]["end"]

    consent = Consent.model_validate(consent_raw)
    assert consent.provision is not None
    assert all(provision.period.end is None for provision in consent.provision.provision)

    research_consent = ResearchConsent(schemaVersion="2026.0.0", scope=consent)
    start = consent.provision.provision[0].period.start.first_moment.date()
    assert ResearchConsent.consents_to_research([research_consent], date=start)
    assert ResearchConsent.consents_to_research([research_consent], date=date(year=2999, month=12, day=31))
    assert not ResearchConsent.consents_to_research([research_consent], date=start - timedelta(days=1)), (
        "a provision must not apply before its start date"
    )


def test_consent_requires_root_provision_period():
    """Every profile version pins the root provision period to 1..1; only its end became optional."""
    consent_raw = _consent_raw("minimal_consented")
    del consent_raw["provision"]["period"]

    with pytest.raises(ValidationError, match="period"):
        Consent.model_validate(consent_raw)


def test_root_provision_period_caps_open_ended_sub_provisions():
    """The root provision frames every nested rule, so a permit must not outlive its period."""
    consent_raw = _consent_raw("minimal_consented")
    for provision in consent_raw["provision"]["provision"]:
        del provision["period"]["end"]

    consent = Consent.model_validate(consent_raw)
    research_consent = ResearchConsent(schemaVersion="2026.0.0", scope=consent)
    root_end_day = consent.provision.period.end.last_moment.date()
    assert ResearchConsent.consents_to_research([research_consent], date=root_end_day), (
        "consent must still apply on the last day of the root provision period"
    )
    assert not ResearchConsent.consents_to_research([research_consent], date=root_end_day + timedelta(days=1)), (
        "an open-ended permit must not outlive the root provision period"
    )


_VERSIONS_REQUIRING_PERIOD_END = tuple(
    version
    for version, profile in RESEARCH_CONSENT_PACKAGE_PROFILES.items()
    if profile in PROFILES_REQUIRING_PERIOD_END
)
# derived from the mapping under test, so an emptied mapping would silently leave nothing to run
assert _VERSIONS_REQUIRING_PERIOD_END, "no accepted package ships a profile requiring a period end"


@pytest.mark.parametrize("version", _VERSIONS_REQUIRING_PERIOD_END)
def test_open_ended_period_rejected_where_the_profile_requires_an_end(version: str):
    """
    A period without an end never expires, so honouring one would outlive what the donor signed.

    The submitter has to be told which periods to fix, not only that one of them is wrong.
    """
    with pytest.raises(ValidationError) as excinfo:
        ResearchConsent(schemaVersion=version, scope=_consent("minimal_consented_open_ended"))

    assert "requires an end on every provision period" in str(excinfo.value)
    assert "provision.period" in str(excinfo.value)
    assert "provision.provision[0].period" in str(excinfo.value)


def test_open_ended_nested_period_alone_is_rejected():
    """The profile pins the end at both levels, so a bounded root does not excuse an open nested one."""
    consent_raw = _consent_raw("minimal_consented")
    for provision in consent_raw["provision"]["provision"]:
        del provision["period"]["end"]

    with pytest.raises(ValidationError, match=r"provision\.provision\[0\]\.period"):
        ResearchConsent(schemaVersion="2025.0.1", scope=Consent.model_validate(consent_raw))


def test_open_ended_period_accepted_without_a_declared_schema_version():
    """Metadata 1.3 onwards may omit the version, which leaves no profile to enforce."""
    research_consent = ResearchConsent(scope=_consent("minimal_consented_open_ended"))
    assert research_consent.scope.provision.period.end is None


def test_research_consent_rejects_category_in_both_spellings():
    """The MII broad consent category must appear exactly once, even across both CodeSystem spellings."""
    consent_raw = _consent_raw("minimal_consented")
    new_style = MII_CONSENT_CATEGORY_SYSTEMS[1]
    mii_category = copy.deepcopy(
        next(c for c in consent_raw["category"] if c["coding"][0]["code"] == MII_BROAD_CONSENT_OID)
    )
    mii_category["coding"][0]["system"] = new_style
    consent_raw["category"].append(mii_category)

    with pytest.raises(ValidationError, match="Duplicate required category"):
        Consent.model_validate(consent_raw)


def test_research_consent_missing_category_message_is_readable():
    """The missing-category error reaches submitters, so it must spell out the accepted system|code pairs."""
    consent_raw = _consent_raw("minimal_consented")
    mii_category = next(c for c in consent_raw["category"] if c["coding"][0]["code"] == MII_BROAD_CONSENT_OID)
    mii_category["coding"][0]["system"] = "https://example.org/not-the-consent-category-system"

    with pytest.raises(ValidationError) as excinfo:
        Consent.model_validate(consent_raw)

    message = str(excinfo.value)
    assert "Missing expected categories" in message
    for system in MII_CONSENT_CATEGORY_SYSTEMS:
        assert f"{system}|{MII_BROAD_CONSENT_OID}" in message
    assert "frozenset" not in message


def test_research_consent_subprovisions_deny_permit():
    """Within one research consent's subprovisions, deny before permit should return a non-consented state."""
    consent_raw = _consent_raw("minimal_nonconsented")

    # add a permit subprovision object for same consent object, after the deny subprovision
    new_permit_subprovision = copy.deepcopy(consent_raw["provision"]["provision"][0])
    new_permit_subprovision["type"] = "permit"
    consent_raw["provision"]["provision"].append(new_permit_subprovision)

    consent = Consent.model_validate(consent_raw)

    assert not ResearchConsent.consents_to_research(
        [ResearchConsent(scope=consent)], date=date(year=2025, month=10, day=13)
    )


def test_research_consents_deny_permit():
    """Having two research consents, where deny comes before permit, should return a non-consented state."""
    consent_raw = _consent_raw("minimal_nonconsented")
    consent1 = Consent.model_validate(consent_raw)

    # add a permit consent object for same donor
    consent_raw["provision"]["provision"][0]["type"] = "permit"
    consent2 = Consent.model_validate(consent_raw)

    assert not ResearchConsent.consents_to_research(
        (ResearchConsent(scope=consent1), ResearchConsent(scope=consent2)), date=date(year=2025, month=10, day=13)
    )


def test_research_consent_no_subprovisions():
    """Consent objects are allowed to have no provisions under the root."""
    consent_json_raw = _consent_raw("minimal_consented")
    del consent_json_raw["provision"]["provision"]
    Consent.model_validate(consent_json_raw)


@pytest.mark.parametrize(
    "case,versions,kinds",
    (
        ("minimal_consented", {BroadConsentVersion.V1_6F}, {ConsentDocumentKind.CONSENT}),
        ("minimal_consented_bc_v1_7_2", {BroadConsentVersion.V1_7_2}, {ConsentDocumentKind.CONSENT}),
        # the MII package itself ships a policy URI without the urn:oid: prefix
        ("minimal_consented_bc_v1_6d_bare_oid", {BroadConsentVersion.V1_6D}, {ConsentDocumentKind.CONSENT}),
        ("withdrawal_complete", {BroadConsentVersion.V1_6F}, {ConsentDocumentKind.COMPLETE_WITHDRAWAL}),
        ("rejection_bc_v1_7_2", {BroadConsentVersion.V1_7_2}, {ConsentDocumentKind.REJECTION}),
    ),
)
def test_consent_derives_broad_consent_version_from_policy_oid(
    case: str, versions: set[BroadConsentVersion], kinds: set[ConsentDocumentKind]
):
    """The signed document version lives in policy[].uri, not in the declared schemaVersion."""
    consent = _consent(case)
    assert consent.broad_consent_versions == versions
    assert consent.document_kinds == kinds
    assert not consent.unknown_document_oids


def test_consent_reports_unknown_policy_oid_without_rejecting_it():
    """A future broad consent version must not break submissions, but must not be silently claimed either."""
    consent = _consent("minimal_consented_unknown_policy_oid")
    assert consent.unknown_document_oids == {"2.16.840.1.113883.3.1937.777.24.2.9999"}
    assert not consent.broad_consent_versions
    assert not consent.document_kinds


def test_consent_derives_version_from_category_coding():
    """
    The version and module OIDs are also valid Consent.category codings, so a consent stating its
    version there rather than only in policy must be recognized too.
    """
    consent_raw = _consent_raw("minimal_consented")
    consent_raw["category"].append(
        {"coding": [{"system": MII_CONSENT_VERSION_MODULES_SYSTEM, "code": "2.16.840.1.113883.3.1937.777.24.2.2079"}]}
    )
    consent = Consent.model_validate(consent_raw)
    assert consent.broad_consent_versions == {BroadConsentVersion.V1_6F, BroadConsentVersion.V1_7_2}


def test_research_consent_exposes_broad_consent_versions():
    """The version is reachable from the submission model without reaching into the FHIR resource."""
    research_consent = ResearchConsent(schemaVersion="2026.0.0", scope=_consent("minimal_consented_bc_v1_7_2"))
    assert research_consent.broad_consent_versions == {BroadConsentVersion.V1_7_2}


def test_research_consent_without_parsed_scope_has_no_broad_consent_version():
    """A scope that fell back to dict must report no version rather than raising."""
    research_consent = ResearchConsent(schemaVersion="2026.0.0", scope={"not": "a consent"})
    assert research_consent.broad_consent_versions == frozenset()


@pytest.mark.parametrize("status", ("draft", "proposed", "rejected", "inactive", "entered-in-error"))
def test_consent_not_in_force_grants_nothing(status: str):
    """
    Consent.status is a FHIR modifier element: anything but 'active' marks the consent as not
    currently valid, so its permit provisions must not be read as research consent.
    """
    consent_raw = _consent_raw("minimal_consented")
    consent_raw["status"] = status
    research_consent = ResearchConsent(schemaVersion="2026.0.0", scope=Consent.model_validate(consent_raw))

    assert not research_consent.scope.is_in_force
    assert research_consent.consent_by_code(date(year=2024, month=1, day=1)) == {}
    assert not ResearchConsent.consents_to_research([research_consent], date=date(year=2024, month=1, day=1))


def test_consent_in_force_still_grants():
    """The status gate must not change the outcome for the active consents it guards."""
    research_consent = ResearchConsent(schemaVersion="2026.0.0", scope=_consent("minimal_consented"))
    assert research_consent.scope.is_in_force
    assert ResearchConsent.consents_to_research([research_consent], date=date(year=2024, month=1, day=1))


def test_consent_not_in_force_warns(caplog):
    """A status that refuses is easy to state by mistake, and otherwise shows up only as a False flag."""
    consent_raw = _consent_raw("minimal_consented")
    consent_raw["status"] = "rejected"

    ResearchConsent(schemaVersion="2026.0.0", scope=Consent.model_validate(consent_raw))

    assert "status 'rejected'" in caplog.text
    assert "grants nothing" in caplog.text


def test_consent_in_force_does_not_warn(caplog):
    """A warning on every valid consent would train submitters to ignore it."""
    ResearchConsent(schemaVersion="2026.0.0", scope=_consent("minimal_consented"))
    assert "grants nothing" not in caplog.text


@pytest.mark.parametrize("case", ("withdrawal_complete", "rejection_bc_v1_7_2", "status_rejected"))
def test_withdrawal_and_rejection_do_not_consent_to_research(case: str):
    research_consent = ResearchConsent(schemaVersion="2026.0.0", scope=_consent(case))
    assert not ResearchConsent.consents_to_research([research_consent], date=date(year=2024, month=1, day=1))


def test_consent_accepts_extra_category_with_several_codings():
    """The category slicing is open, so a category outside the pinned slices may carry several codings."""
    consent = _consent("minimal_consented_extra_multi_coding_category")
    assert len(consent.category) == 3
    assert len(consent.category[2].coding) == 2


def test_consent_rejects_several_codings_on_the_mii_category():
    """The mii slice itself is pinned to exactly one coding."""
    with pytest.raises(ValidationError, match="carries the mii category"):
        _consent("invalid_mii_category_multi_coding")


def test_consent_rejects_provisions_it_would_not_evaluate():
    """Both carry permissions consent evaluation never reads, so accepting them would lose them."""
    with pytest.raises(ValidationError, match=re.escape("provision.provision[].provision is not allowed")):
        _consent("invalid_third_level_provision")
    with pytest.raises(ValidationError, match=re.escape("consent.provision.code is not allowed")):
        _consent("invalid_root_provision_code")


def test_consent_rejects_other_resource_types():
    """Unknown fields are ignored, so without this guard any object with consent-shaped keys parses."""
    with pytest.raises(ValidationError, match="Input should be 'Consent'"):
        _consent("invalid_wrong_resource_type")


def test_consent_keeps_patient_identifier():
    """The profile allows identifying the patient by identifier; dropping it would lose the pseudonym."""
    consent = _consent("minimal_consented_patient_by_identifier")
    assert consent.patient.reference is None
    assert consent.patient.identifier is not None
    assert consent.patient.identifier.value == "42"


def test_consent_rejects_unidentified_patient():
    """
    The profile requires neither reference nor identifier (both are mustSupport, 0..1), but a
    patient stating neither, e.g. by display name only, carries nothing this model keeps.
    """
    consent_raw = _consent_raw("minimal_consented")
    consent_raw["patient"] = {}

    with pytest.raises(ValidationError, match="reference or an identifier"):
        Consent.model_validate(consent_raw)


def test_date_only_provision_end_covers_the_whole_day():
    """FHIR treats a date-only end as inclusive of that day, not as midnight."""
    consent = _consent("minimal_consented")
    provision = consent.provision.provision[0]
    end_day = provision.period.end.last_moment.date()

    assert provision.period.contains(datetime(end_day.year, end_day.month, end_day.day, 12, 0, tzinfo=UTC))
    assert provision.period.contains(datetime(end_day.year, end_day.month, end_day.day, 23, 59, tzinfo=UTC))
    assert not provision.period.contains(
        datetime(end_day.year, end_day.month, end_day.day, 0, 0, tzinfo=UTC) + timedelta(days=1)
    )

    research_consent = ResearchConsent(schemaVersion="2026.0.0", scope=consent)
    noon = datetime(end_day.year, end_day.month, end_day.day, 12, 0, tzinfo=UTC)
    assert research_consent.consent_by_code(noon), "consent must still apply at noon on its last valid day"


def test_datetime_provision_end_is_not_extended_to_the_whole_day():
    """A period end that states a time must keep it, otherwise it silently gains up to a day."""
    consent = _consent("minimal_consented_with_datetime")
    berlin_summer = timezone(timedelta(hours=2))
    assert consent.provision.period.end.last_moment == datetime(2050, 8, 31, 18, 7, tzinfo=berlin_summer)
    period = consent.provision.provision[0].period
    assert period.end.last_moment == datetime(2025, 8, 31, 9, 4, 51, tzinfo=berlin_summer)
    assert not period.contains(datetime(2025, 8, 31, 12, 0, tzinfo=UTC))


@pytest.mark.parametrize(
    ("bound", "start", "end"),
    [
        ("2020", datetime(2020, 1, 1, 0, 0, tzinfo=UTC), datetime(2020, 12, 31, 23, 59, 59, 999999, tzinfo=UTC)),
        ("2020-09", datetime(2020, 9, 1, 0, 0, tzinfo=UTC), datetime(2020, 9, 30, 23, 59, 59, 999999, tzinfo=UTC)),
        ("2020-02", datetime(2020, 2, 1, 0, 0, tzinfo=UTC), datetime(2020, 2, 29, 23, 59, 59, 999999, tzinfo=UTC)),
        ("2021-02", datetime(2021, 2, 1, 0, 0, tzinfo=UTC), datetime(2021, 2, 28, 23, 59, 59, 999999, tzinfo=UTC)),
        ("2020-09-01", datetime(2020, 9, 1, 0, 0, tzinfo=UTC), datetime(2020, 9, 1, 23, 59, 59, 999999, tzinfo=UTC)),
        # the year regex spells out its bounds by hand, so exercise both ends of the range it allows
        ("0001", datetime(1, 1, 1, 0, 0, tzinfo=UTC), datetime(1, 12, 31, 23, 59, 59, 999999, tzinfo=UTC)),
        ("9999", datetime(9999, 1, 1, 0, 0, tzinfo=UTC), datetime(9999, 12, 31, 23, 59, 59, 999999, tzinfo=UTC)),
    ],
)
def test_period_bound_covers_the_whole_span_it_names(bound: str, start: datetime, end: datetime):
    """
    FHIR dateTime states a year or a month in place of a full date, each naming a whole span.

    A start therefore begins at that span's first moment and an end expires at its last, the same
    rule FHIR spells out for a date-only end.
    """
    assert Period(start=bound).start.first_moment == start
    assert Period(start="2020-01-01", end=bound).end.last_moment == end


def test_period_bound_of_reduced_precision_stays_in_force_for_its_whole_span():
    """
    A permit ending in a stated year must run to that year's end.

    Read as a Unix timestamp instead, such a bound lands in 1970 and refuses consent it grants.
    """
    consent_raw = _consent_raw("minimal_consented")
    consent_raw["provision"]["period"]["end"] = "2030"
    for provision in consent_raw["provision"]["provision"]:
        provision["period"]["end"] = "2030"

    research_consent = ResearchConsent(schemaVersion="2026.0.0", scope=Consent.model_validate(consent_raw))
    assert ResearchConsent.consents_to_research([research_consent], date=date(year=2024, month=1, day=1))
    assert ResearchConsent.consents_to_research([research_consent], date=date(year=2030, month=12, day=31))
    assert not ResearchConsent.consents_to_research([research_consent], date=date(year=2031, month=1, day=1))


@pytest.mark.parametrize(
    ("bound", "expected"),
    [
        ("2020-09-01T14:37:22Z", datetime(2020, 9, 1, 14, 37, 22, tzinfo=UTC)),
        (
            "2020-09-01T14:37:22.123456+02:00",
            datetime(2020, 9, 1, 14, 37, 22, 123456, tzinfo=timezone(timedelta(hours=2))),
        ),
        ("2020-09-01T14:37:22-05:00", datetime(2020, 9, 1, 14, 37, 22, tzinfo=timezone(timedelta(hours=-5)))),
        # +14:00 and -13:59 are the outermost offsets the regex allows, each spelled out separately
        ("2020-09-01T00:00:00+14:00", datetime(2020, 9, 1, 0, 0, tzinfo=timezone(timedelta(hours=14)))),
        (
            "2020-09-01T23:59:59-13:59",
            datetime(2020, 9, 1, 23, 59, 59, tzinfo=timezone(-timedelta(hours=13, minutes=59))),
        ),
    ],
)
def test_period_bound_carrying_a_time_is_kept_as_submitted(bound: str, expected: datetime):
    """A bound that already states a time names one moment, so neither end of a period widens it."""
    start = Period(start=bound).start
    end = Period(start="2020-01-01", end=bound).end

    assert str(start) == bound, "the submitted spelling must survive"
    assert start.first_moment == expected
    assert end.last_moment == expected


@pytest.mark.parametrize("bound", ["2020-09-01T00:00:00", "2020-09-01T14:37:22", "2020-12-31T23:59:59.999999"])
def test_period_bound_without_a_timezone_is_read_but_not_conformant(bound: str):
    """
    A bound written before FHIR's timezone rule was enforced still parses, so an old document
    stays readable, and it still resolves to a moment by reading the zone as UTC.

    Refusing it is a submission-time rule, checked where a submission is validated, so nothing
    here has to guess on the submitter's behalf or rewrite what they wrote.
    """
    start = Period(start=bound).start

    assert str(start) == bound
    assert not start.is_valid_fhir
    assert start.first_moment.tzinfo is not None


@pytest.mark.parametrize(
    ("value", "states_a_time", "has_timezone", "is_valid_fhir"),
    [
        ("2020", False, False, True),  # no time, so no zone is required
        ("2020-09", False, False, True),
        ("2020-09-01", False, False, True),
        ("2020-09-01T14:37:22Z", True, True, True),
        ("2020-09-01T14:37:22+02:00", True, True, True),
        ("2020-09-01T14:37:22", True, False, False),  # the one shape FHIR forbids
    ],
)
def test_fhir_datetime_reports_what_it_names(value: str, states_a_time, has_timezone, is_valid_fhir):
    """
    The pattern admits a superset of FHIR, so a parsed value has to say whether FHIR permits it.

    Shape and validity are separate questions: a zone-less time parses, which is what allows it to
    be repaired, but it is not valid FHIR.
    """
    parsed = FhirDateTime.parse(value)

    assert parsed.raw == value
    assert str(parsed) == value
    assert parsed.states_a_time is states_a_time
    assert parsed.has_timezone is has_timezone
    assert parsed.is_valid_fhir is is_valid_fhir


@pytest.mark.parametrize("mode", ["python", "json"])
def test_period_dumps_a_bound_as_the_value_it_is(mode: str, recwarn):
    """
    Dumping must hand back the value in either mode, not this type's own fields.

    A python-mode dump that serialized the dataclass would leak `raw`, `year` and the rest to any
    caller not passing mode="json", and warn while doing it.
    """
    dumped = Period(start="2020-09-01").model_dump(mode=mode)

    assert dumped["start"] == "2020-09-01"
    assert not [warning for warning in recwarn if "Serialization" in str(warning.message)]


@pytest.mark.parametrize("value", ["20200901", "2020-09-01 14:37:22", "1700000000", "", 1700000000, None])
def test_fhir_datetime_refuses_what_is_not_shaped_like_one(value):
    """Anything outside the pattern has no reading at all, so it does not parse."""
    with pytest.raises(ValueError, match="is not a FHIR dateTime"):
        FhirDateTime.parse(value)


@pytest.mark.parametrize(
    ("value", "first", "last"),
    [
        ("2020", date(2020, 1, 1), date(2020, 12, 31)),
        ("2020-02", date(2020, 2, 1), date(2020, 2, 29)),
        ("2021-02", date(2021, 2, 1), date(2021, 2, 28)),
        ("2020-09-01", date(2020, 9, 1), date(2020, 9, 1)),
    ],
)
def test_fhir_datetime_covers_the_span_its_precision_names(value: str, first: date, last: date):
    parsed = FhirDateTime.parse(value)

    assert parsed.covered_days() == (first, last)
    assert parsed.first_moment == datetime.combine(first, datetime.min.time(), tzinfo=UTC)
    assert parsed.last_moment.date() == last


@pytest.mark.parametrize("value", ["2021-02-30", "0000", "2020-04-31"])
def test_fhir_datetime_refuses_a_date_that_does_not_exist(value: str):
    """
    The pattern admits any day up to 31 in any month, so parsing checks the date exists.

    Catching it here means the document is still in hand to blame, rather than something later
    failing when it asks the value for a moment.
    """
    assert FhirDateTime._PATTERN.fullmatch(value), "the shape is fine; only the date is not"

    with pytest.raises(ValueError, match="is not an existing date"):
        FhirDateTime.parse(value)


def test_period_bound_with_a_leap_second_cannot_be_resolved():
    """
    Known gap: FHIR admits a leap second and the pattern follows it, but datetime cannot hold one.

    Such a bound parses, so a document carrying one is still readable, and only resolving it to a
    moment fails.
    """
    start = Period(start="2020-09-01T23:59:60Z").start

    assert str(start) == "2020-09-01T23:59:60Z"
    with pytest.raises(ValueError, match="cannot be represented"):
        _ = start.first_moment


def test_period_bound_given_as_a_python_object_is_spelled_as_a_document_would():
    """
    A model built in code hands the field a date or datetime rather than the string a document has.

    Those are spelled the way a document would state them, so a model built in code and one parsed
    from JSON hold the same value.
    """
    assert str(Period(start=date(2020, 9, 1)).start) == "2020-09-01"
    assert Period(start=date(2020, 9, 1)).start.first_moment == datetime(2020, 9, 1, 0, 0, tzinfo=UTC)
    assert Period(start=date(2020, 1, 1), end=date(2020, 9, 1)).end.last_moment == datetime(
        2020, 9, 1, 23, 59, 59, 999999, tzinfo=UTC
    )

    berlin_summer = timezone(timedelta(hours=2))
    aware = Period(start=datetime(2020, 9, 1, 14, 37, 0, tzinfo=berlin_summer)).start
    assert str(aware) == "2020-09-01T14:37:00+02:00"
    assert aware.first_moment == datetime(2020, 9, 1, 12, 37, tzinfo=UTC)


def test_period_bound_refuses_a_naive_datetime():
    """
    A datetime built in code without a zone is refused rather than read as UTC.

    Guessing is what this type exists to avoid, and code building a consent is the one place where
    the caller still knows which zone was meant. Reading it as UTC would instead produce a value
    that only fails later, when the submission is validated.
    """
    with pytest.raises(ValidationError, match="names no timezone"):
        Period(start=datetime(2020, 9, 1, 14, 37))


def test_period_bound_accepts_a_plain_date():
    """A date names no time, so FHIR asks for no zone and nothing has to be guessed."""
    assert str(Period(start=date(2020, 9, 1)).start) == "2020-09-01"


@pytest.mark.parametrize(
    "bound",
    [
        "2020-09-01 14:37:22+02:00",  # FHIR separates date and time by 'T', never by a space
        "2020-09-01T14:37",  # seconds are not optional, unlike the timezone
        "20200901",  # the basic format is not a FHIR dateTime
        "2020-W36-2",  # nor is a week date
        "2020-9",  # a month is two digits
        "20",  # a year is four digits, so a truncated one must not read as year 20
        "202",
        "20200",  # nor may a fifth digit be ignored
        "2020-00",  # month zero does not exist
        "2020-13",  # nor does month 13
        "2020-09-00",  # nor day zero
        "2020-09-32",  # nor day 32
        "2020-09-01T24:00:00Z",  # hours run to 23; FHIR has no end-of-day 24:00
        "2020-09-01T14:60:00Z",  # minutes run to 59, unlike seconds, which allow a leap second
        "2020-09-01T14:37:22+15:00",  # no timezone is that far ahead
        "2020-09-01T14:37:22+14:01",  # +14:00 is the furthest offset, and only exactly on the hour
        "2020-09-01T14:37:22+2:00",  # an offset hour is two digits
        "1700000000",  # a Unix timestamp would otherwise be read as one
        1700000000,
        "",  # an empty string names nothing
        "2020-09-01\n",  # a trailing newline must not slip past the anchoring
    ],
)
def test_period_rejects_a_bound_that_is_not_a_fhir_datetime(bound):
    """
    Pydantic reads an all-numeric string as a Unix timestamp, placing such a bound in 1970 unremarked
    and refusing consent that was granted, so anything outside the dateTime regex must be rejected.

    The match is exact rather than an alternation over both rejection messages, so that a change
    rerouting a value to the other message shows up here instead of passing quietly.
    """
    with pytest.raises(ValidationError, match="is not a FHIR dateTime"):
        Period(start=bound)
    with pytest.raises(ValidationError, match="is not a FHIR dateTime"):
        Period(start="2020-01-01", end=bound)


@pytest.mark.parametrize(
    "bound",
    [
        "2021-02-30",  # the regex admits any day up to 31, but February has no 30th
        "2020-02-30",  # not even in a leap year
        "2020-04-31",  # April has 30 days
        "0000",  # the regex admits any four digits, but there is no year zero
        "0000-01-01",
    ],
)
def test_period_rejects_a_bound_that_names_no_existing_date(bound: str):
    """
    The regex checks shape, not whether the date exists, so `_covered_days` rejects the rest.

    These reach a different error than a malformed bound does, and the split keeps that distinction
    visible: were one path to swallow the other's cases, one of these two tests would fail.
    """
    with pytest.raises(ValidationError, match="is not an existing date"):
        Period(start=bound)
    with pytest.raises(ValidationError, match="is not an existing date"):
        Period(start="2020-01-01", end=bound)


def test_period_bound_is_revalidated_on_assignment():
    """
    The base model sets validate_assignment, so a bound replaced after construction is checked too.

    Without this the widening and the timestamp rejection would apply only at construction, and code
    that patches a period in place could reintroduce the 1970 bug the regex exists to prevent.
    """
    period = Period(start="2020-01-01")

    period.end = "2030"
    assert period.end.last_moment == datetime(2030, 12, 31, 23, 59, 59, 999999, tzinfo=UTC)

    with pytest.raises(ValidationError, match="is not a FHIR dateTime"):
        period.end = "1700000000"


@pytest.mark.parametrize("case", VALID_CONSENT_CASES)
def test_consent_round_trips_through_its_own_json(case: str):
    """
    A serialized consent must still be a document the model accepts.

    FHIR demands a timezone alongside a time, so a bound widened to a datetime has to carry one for
    its serialization to stay a dateTime.
    """
    consent = _consent(case)
    assert Consent.model_validate_json(consent.model_dump_json(by_alias=True)) == consent


@pytest.mark.parametrize("case", VALID_CONSENT_CASES)
def test_consent_stored_as_naive_datetimes_still_parses(case: str):
    """
    Stored submission_metadata was dumped from the naive datetimes the previous validator produced,
    so every date-only value in the database reads "2020-09-01T00:00:00".

    Keeping the submitted value means such a document parses as it stands: `db backfill` and every
    other reader of an archived submission need no repair step, and nothing has to be migrated.
    """
    stored = _as_stored_by_the_previous_validator(_consent_raw(case))

    consent = Consent.model_validate(stored)

    assert consent.date_time.first_moment.tzinfo is not None
    assert Consent.model_validate_json(consent.model_dump_json(by_alias=True)) == consent


def _as_stored_by_the_previous_validator(raw):
    """Rewrite every date-only value the way a naive `model_dump(mode="json")` would have."""
    if isinstance(raw, dict):
        return {key: _as_stored_by_the_previous_validator(value) for key, value in raw.items()}
    if isinstance(raw, list):
        return [_as_stored_by_the_previous_validator(value) for value in raw]
    if isinstance(raw, str) and re.fullmatch(r"[0-9]{4}-[0-9]{2}-[0-9]{2}", raw):
        return f"{raw}T00:00:00"
    return raw


@pytest.mark.parametrize("case", VALID_CONSENT_CASES)
def test_consent_round_trips_without_rewriting_what_was_submitted(case: str):
    """
    Serializing must return the document, not a normalised reading of it.

    A date that came back as a midnight timestamp would lose the fact that it named a whole day,
    which is exactly how stored metadata came to hold ends that expire a day early.
    """
    raw = _consent_raw(case)

    dumped = json.loads(Consent.model_validate(raw).model_dump_json(by_alias=True, exclude_none=True))

    assert dumped["dateTime"] == raw["dateTime"]
    if period := raw.get("provision", {}).get("period"):
        assert dumped["provision"]["period"]["start"] == period["start"]


def test_consent_reports_the_datetimes_fhir_does_not_permit():
    """
    The report walks the parsed model, so a period buried in a nested provision cannot hide one and
    no raw document has to be re-read to find it.
    """
    raw = _consent_raw("minimal_consented")
    raw["dateTime"] = "2020-09-01T14:37:22"
    raw["provision"]["period"]["end"] = "2050-08-31T00:00:00"
    raw["provision"]["provision"][0]["period"]["start"] = "2020-09-01T08:00:00+02:00"
    raw["provision"]["provision"][0]["period"]["end"] = "2025-01-01T00:00:00"
    raw["verification"] = [{"verified": True, "verificationDate": "2021-01-01T09:00:00"}]

    offending = [str(value) for value in Consent.model_validate(raw).datetimes_fhir_does_not_permit()]

    assert offending == [
        "2020-09-01T14:37:22",
        "2050-08-31T00:00:00",
        "2025-01-01T00:00:00",
        "2021-01-01T09:00:00",
    ]


@pytest.mark.parametrize("case", VALID_CONSENT_CASES)
def test_shipped_consent_examples_are_fhir_conformant(case: str):
    """The bundled examples must stay conformant, or they would fail validation on submission."""
    assert Consent.model_validate(_consent_raw(case)).datetimes_fhir_does_not_permit() == []


def test_consent_date_time_rejects_a_value_that_is_not_a_fhir_datetime():
    """Consent.dateTime is a FHIR dateTime like the period bounds and takes the same timestamp path."""
    consent_raw = _consent_raw("minimal_consented")
    consent_raw["dateTime"] = "1700000000"

    with pytest.raises(ValidationError, match="FHIR dateTime"):
        Consent.model_validate(consent_raw)


def test_consent_date_time_of_reduced_precision_takes_its_first_moment():
    """A point in time bounds nothing, so a year names the first moment it can denote."""
    consent_raw = _consent_raw("minimal_consented")
    consent_raw["dateTime"] = "2020"

    assert Consent.model_validate(consent_raw).date_time.first_moment == datetime(2020, 1, 1, 0, 0, tzinfo=UTC)


def test_consent_date_time_carrying_a_time_is_kept_as_submitted():
    """Widening applies only to a value without a time; a full dateTime must survive its validator."""
    consent_raw = _consent_raw("minimal_consented")
    consent_raw["dateTime"] = "2020-09-01T14:37:22+02:00"

    assert Consent.model_validate(consent_raw).date_time.first_moment == datetime(
        2020, 9, 1, 14, 37, 22, tzinfo=timezone(timedelta(hours=2))
    )


def test_verification_date_rejects_a_value_that_is_not_a_fhir_datetime():
    """Consent.verification.verificationDate is a FHIR dateTime and takes the same timestamp path."""
    consent_raw = _consent_raw("minimal_consented")
    consent_raw["verification"] = [{"verified": True, "verificationDate": "1700000000"}]

    with pytest.raises(ValidationError, match="FHIR dateTime"):
        Consent.model_validate(consent_raw)


def test_verification_date_of_reduced_precision_takes_its_first_moment():
    """Like Consent.dateTime, a verification date is a point in time and bounds nothing."""
    consent_raw = _consent_raw("minimal_consented")
    consent_raw["verification"] = [{"verified": True, "verificationDate": "2020"}]

    assert Consent.model_validate(consent_raw).verification[0].verification_date.first_moment == datetime(
        2020, 1, 1, 0, 0, tzinfo=UTC
    )


def test_verification_date_carrying_a_time_is_kept_as_submitted():
    """Widening applies only to a value without a time; a full dateTime must survive its validator."""
    consent_raw = _consent_raw("minimal_consented")
    consent_raw["verification"] = [{"verified": True, "verificationDate": "2020-09-01T14:37:22+02:00"}]

    assert Consent.model_validate(consent_raw).verification[0].verification_date.first_moment == datetime(
        2020, 9, 1, 14, 37, 22, tzinfo=timezone(timedelta(hours=2))
    )


def test_verification_date_stays_optional():
    """FHIR marks verificationDate as 0..1, so a verification without one must still validate."""
    consent_raw = _consent_raw("minimal_consented")
    consent_raw["verification"] = [{"verified": True}]

    assert Consent.model_validate(consent_raw).verification[0].verification_date is None


### Conformance against the MII artefacts vendored in example_terminology. These tests discover
### their inputs, so dropping a newly released version in adds cases. See that directory's README.


def _vendored(prefix: str) -> list[str]:
    """File names of every vendored MII artefact whose name starts with ``prefix``."""
    names = sorted(
        entry.name
        for entry in importlib.resources.files(example_terminology).iterdir()
        if entry.name.startswith(prefix) and entry.name.endswith(".json")
    )
    assert names, f"no vendored MII artefact matches {prefix!r}"
    return names


def _load_vendored(name: str) -> dict:
    return json.loads(importlib.resources.files(example_terminology).joinpath(name).read_text())


def _package_manifest() -> dict[str, dict[str, str]]:
    """Which version each artefact of a published package states, as the vendoring script recorded it."""
    return json.loads(importlib.resources.files(example_terminology).joinpath("packages.json").read_text())


def _declared_version(name: str) -> str:
    """
    The version a vendored file name claims, e.g. 'mii-cs-consent-policy.1.1.0.json' -> '1.1.0'.

    Artefacts are named after their own resource version rather than the package that ships them:
    several packages ship the same CodeSystem or profile version, and the mapping between the two
    is recorded in packages.json.
    """
    return name.removesuffix(".json").split(".", 1)[1]


def _artefact_prefix(name: str) -> str:
    """The artefact a vendored file name claims, e.g. 'mii-cs-consent-policy.1.1.0.json' -> 'mii-cs-consent-policy'."""
    return name.removesuffix(".json").split(".", 1)[0]


def _codesystem_concepts(codesystem: dict) -> dict[str, dict[str, object]]:
    """Every concept of a CodeSystem, including nested ones, mapped to its properties."""
    concepts: dict[str, dict[str, object]] = {}

    def collect(nodes):
        for node in nodes:
            concepts[node["code"]] = {
                prop["code"]: prop.get("valueString", prop.get("valueBoolean", prop.get("valueCode")))
                for prop in node.get("property", [])
            }
            collect(node.get("concept", []))

    collect(codesystem.get("concept", []))
    return concepts


@pytest.mark.parametrize("name", _vendored("mii-"))
def test_vendored_artefact_holds_the_version_its_name_claims(name: str):
    """The file name is what the other tests report, so it must not drift from the contents."""
    assert _load_vendored(name)["version"] == _declared_version(name)


@pytest.mark.parametrize("name", _vendored("mii-cs-consent-version-modules."))
def test_document_oid_table_covers_every_shipped_oid(name: str):
    """Every OID the MII ships must be classified, so a new broad consent version fails here."""
    shipped = {concept["code"] for concept in _load_vendored(name)["concept"]}
    missing = shipped - set(BROAD_CONSENT_DOCUMENT_OIDS)
    assert not missing, f"{name} defines OIDs the model does not classify: {sorted(missing)}"


def test_document_oid_table_invents_no_oids():
    """Every OID claimed by the table must come from some vendored CodeSystem, not from guesswork."""
    shipped = {
        concept["code"]
        for name in _vendored("mii-cs-consent-version-modules.")
        for concept in _load_vendored(name)["concept"]
    }
    assert set(BROAD_CONSENT_DOCUMENT_OIDS) <= shipped


@pytest.mark.parametrize("name", _vendored("mii-cs-consent-policy."))
def test_evaluated_policy_codes_are_active(name: str):
    """
    The two codes research consent is derived from must exist and stay active in every vendored
    policy CodeSystem, otherwise consent would be granted on a permission the MII has withdrawn.
    """
    concepts = _codesystem_concepts(_load_vendored(name))
    for code in ResearchConsentCodes:
        properties = concepts.get(code.value)
        assert properties is not None, f"{code.value} is missing from {name}"
        assert properties.get("status") != "deprecated", f"{code.value} is deprecated in {name}"
        assert not properties.get("inactive"), f"{code.value} is inactive in {name}"


def _concept_subtree(nodes: list[dict], code: str) -> dict | None:
    """The concept carrying ``code``, searched at any depth."""
    for node in nodes:
        if node["code"] == code:
            return node
        if (found := _concept_subtree(node.get("concept", []), code)) is not None:
            return found
    return None


@pytest.mark.parametrize("name", _vendored("mii-cs-consent-policy."))
def test_scientific_use_stays_within_the_patdat_module(name: str):
    """
    consents_to_research accepts a permit of either research code alone. That is sound only while
    the MII keeps the scientific-use permission a descendant of the PATDAT module, so that a permit
    of the module subsumes it.
    """
    module = _concept_subtree(_load_vendored(name)["concept"], ResearchConsentCodes.PATDAT_ERHEBEN_SPEICHERN_NUTZEN)
    assert module is not None, f"{name} does not define the PATDAT module"
    assert ResearchConsentCodes.MDAT_WISSENSCHAFTLICH_NUTZEN_EU_DSGVO_NIVEAU in _codesystem_concepts(module)


def _differential(profile: dict) -> dict[str, dict]:
    return {element["id"]: element for element in profile["differential"]["element"]}


def _category_slice_coding(elements: dict[str, dict], slice_name: str) -> tuple[str, str]:
    """
    The (system, code) a category slice is pinned to.

    Profile 1.0.8 splits the pattern across the slice and its coding; 1.0.9 states both on the slice.
    """
    coding: dict[str, str] = {}
    for element_id in (f"Consent.category:{slice_name}", f"Consent.category:{slice_name}.coding"):
        element = elements.get(element_id, {})
        pattern = element.get("patternCodeableConcept", {}).get("coding", [{}])[0] | element.get("patternCoding", {})
        coding |= {key: value for key, value in pattern.items() if key in ("system", "code")}
    return coding["system"], coding["code"]


@pytest.mark.parametrize("name", _vendored("mii-pr-consent-einwilligung."))
def test_model_matches_the_profile(name: str):
    """What the model hard-codes about the profile, checked against the profile itself."""
    elements = _differential(_load_vendored(name))

    assert elements["Consent.scope.coding.system"]["fixedUri"] == EXPECTED_SCOPE_CODING_SYSTEM
    assert elements["Consent.scope.coding.code"]["fixedCode"] == EXPECTED_SCOPE_CODING_CODE

    assert _category_slice_coding(elements, "loinc") in EXPECTED_CATEGORIES["loinc"]
    assert _category_slice_coding(elements, "mii") in EXPECTED_CATEGORIES["mii"]
    assert elements["Consent.category"]["min"] == Consent.model_fields["category"].metadata[0].min_length

    # the model rejects what the profile forbids, because evaluation would otherwise drop it silently
    assert elements["Consent.provision.code"]["max"] == "0"
    assert elements["Consent.provision.provision.provision"]["max"] == "0"

    # the model must demand a period end exactly where the profile does, at both provision levels
    requires_end = _declared_version(name) in PROFILES_REQUIRING_PERIOD_END
    for element_id in ("Consent.provision.period.end", "Consent.provision.provision.period.end"):
        assert (elements[element_id]["min"] == 1) == requires_end, (
            f"{name} and PROFILES_REQUIRING_PERIOD_END disagree on {element_id}"
        )
    # the field itself stays optional, since one model serves every profile version
    assert not Period.model_fields["end"].is_required()

    for element_id in ("Consent.provision.period.start", "Consent.provision.provision.period.start"):
        assert elements[element_id]["min"] == 1
    assert Period.model_fields["start"].is_required()

    # the period itself stays required at both provision levels, framing consent evaluation
    for element_id in ("Consent.provision.period", "Consent.provision.provision.period"):
        assert elements[element_id]["min"] == 1
    assert RootConsentProvision.model_fields["period"].is_required()
    assert ConsentProvision.model_fields["period"].is_required()

    # a patient identifier, when given, must carry system and value
    for element_id in ("Consent.patient.identifier.system", "Consent.patient.identifier.value"):
        assert elements[element_id]["min"] == 1
    assert Identifier.model_fields["system"].is_required()
    assert Identifier.model_fields["value"].is_required()


def test_recorded_and_vendored_artefacts_agree():
    """
    A vendored file no package lists would have the model checked against a version no accepted
    schemaVersion can select; a version a package lists but nothing vendors would leave that
    package checked against no artefact at all.
    """
    vendored = {(_artefact_prefix(name), _declared_version(name)) for name in _vendored("mii-")}
    recorded = {
        (prefix, version) for artefacts in _package_manifest().values() for prefix, version in artefacts.items()
    }

    assert recorded == vendored


def test_accepted_packages_ship_the_profiles_the_model_assumes():
    """
    The profile a package ships decides which cardinalities apply, and it is stated by the package
    rather than by the artefacts, which name only themselves. The mapping must repeat what the
    vendoring script recorded, so that a newly released package cannot be classified from memory.
    """
    shipped = {package: artefacts["mii-pr-consent-einwilligung"] for package, artefacts in _package_manifest().items()}
    assert RESEARCH_CONSENT_PACKAGE_PROFILES == shipped


def test_research_consent_schema_version_json_schema_states_the_accepted_versions():
    """
    The published schema constrains schemaVersion only through the restated enum, so a dropped
    restatement leaves it an unconstrained string.
    """
    schema = GrzSubmissionMetadata.model_json_schema()
    # the field is optional, so pydantic wraps the declared schema in an anyOf with null
    schema_version = schema["$defs"]["ResearchConsent"]["properties"]["schemaVersion"]
    enums = [branch["enum"] for branch in schema_version["anyOf"] if "enum" in branch]
    assert enums == [list(RESEARCH_CONSENT_SCHEMA_VERSIONS)]


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
    metadata = json.loads(_metadata_raw(dataset, version))
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
    metadata = json.loads(_metadata_raw(dataset, version))
    metadata["submission"]["diseaseType"] = DiseaseType.rare.value
    metadata["submission"]["submissionDate"] = "2026-05-31"

    GrzSubmissionMetadata.model_validate_json(json.dumps(metadata))

    assert "starting 01.06.2026" in caplog.text


@pytest.mark.parametrize("version", TESTED_VERSIONS)
def test_disease_type_rare_index_has_wgs_and_panel_passes(version: str):
    """
    If the index donor has at least WGS but additionally some panel data, it should still pass.
    """
    metadata = json.loads(_metadata_raw("wgs_trio", version))
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

    GrzSubmissionMetadata.model_validate_json(json.dumps(metadata))


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
    metadata = json.loads(_metadata_raw(dataset, version))
    metadata["submission"]["diseaseType"] = DiseaseType.rare.value
    metadata["submission"]["submissionDate"] = "2026-06-01"

    GrzSubmissionMetadata.model_validate_json(json.dumps(metadata))


@pytest.mark.parametrize("version", TESTED_VERSIONS)
def test_disease_type_rare_non_index_unrestricted(version: str):
    """Non-index donors (e.g., parents in a trio) are exempt from the WGS rule."""
    metadata = json.loads(_metadata_raw("wgs_trio", version))
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

    # should pass because the index donor still has WGS
    GrzSubmissionMetadata.model_validate_json(json.dumps(metadata))


@pytest.mark.parametrize(
    "version,justification",
    itertools.product(
        TESTED_VERSIONS,
        [ResearchConsentNoScopeJustification.LE_TECH.value, ResearchConsentNoScopeJustification.LE_ORG.value],
    ),
)
def test_no_scope_justification_tech_org_raises(version: str, justification: str):
    """LE_TECH and LE_ORG justifications should raise an error starting exactly 01.06.2026."""
    metadata = json.loads(_metadata_raw("wgs_trio", version))
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
    metadata = json.loads(_metadata_raw("wgs_trio", version))
    metadata["submission"]["submissionDate"] = "2026-05-31"

    for donor in metadata["donors"]:
        for consent in donor.get("researchConsents", []):
            consent["noScopeJustification"] = justification
            consent.pop("scope", None)

    GrzSubmissionMetadata.model_validate_json(json.dumps(metadata))

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
    metadata = json.loads(_metadata_raw("wgs_trio", version))
    metadata["submission"]["submissionDate"] = "2026-06-01"

    for donor in metadata["donors"]:
        for consent in donor.get("researchConsents", []):
            consent["noScopeJustification"] = justification
            consent.pop("scope", None)

    GrzSubmissionMetadata.model_validate_json(json.dumps(metadata))
