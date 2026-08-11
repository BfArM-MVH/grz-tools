from datetime import date, datetime, time
from enum import StrEnum
from typing import Annotated, Any, Literal

from pydantic import ConfigDict, Field, field_validator, model_validator

from ..common import StrictBaseModel, as_aware_datetime


class StrictIgnoringBaseModel(StrictBaseModel):
    model_config = ConfigDict(extra="ignore")


class Extension(StrictBaseModel):
    # Extensions carry an open-ended value[x] payload, so extra keys are allowed here only.
    model_config = ConfigDict(extra="allow")
    url: str


class FhirElement(StrictBaseModel):
    # Every FHIR element may carry an id and extensions; any other unknown field stays forbidden.
    id: str | None = None
    extension: list[Extension] | None = None


class ProvisionType(StrEnum):
    DENY = "deny"
    PERMIT = "permit"


class Status(StrEnum):
    DRAFT = "draft"
    PROPOSED = "proposed"
    ACTIVE = "active"
    REJECTED = "rejected"
    INACTIVATE = "inactive"
    ENTERED_IN_ERROR = "entered-in-error"


def _date_to_datetime(value: Any, time_of_day: time) -> Any:
    """
    Widen a date-only period bound to a datetime at ``time_of_day``.

    Values that already carry a time component are passed through untouched for pydantic to parse.

    :param value: raw period bound as submitted.
    :param time_of_day: time to combine a date-only bound with.
    :returns: a datetime for date-only bounds, otherwise ``value`` unchanged.
    """
    if isinstance(value, datetime):
        return value

    if isinstance(value, date):
        return datetime.combine(value, time_of_day)

    if isinstance(value, str):
        has_time_component = "T" in value or " " in value
        if has_time_component:
            return value
        try:
            parsed = date.fromisoformat(value)
        except ValueError:
            return value
        return datetime.combine(parsed, time_of_day)

    return value


class Period(FhirElement):
    start: datetime
    #: Optional because one model serves every profile version and 1.0.9 relaxed this to 0..1. FHIR
    #: reads a period without an end as still running, so ResearchConsent rejects one when the
    #: declared schemaVersion ships a profile that still pins the end to 1..1.
    end: datetime | None = None

    # FHIR reads a date-only bound as the whole day: a start begins at midnight, an end expires
    # at the end of that day.

    @field_validator("start", mode="before")
    @classmethod
    def start_of_day(cls, v):
        return _date_to_datetime(v, time.min)

    @field_validator("end", mode="before")
    @classmethod
    def end_of_day(cls, v):
        return _date_to_datetime(v, time.max)

    def contains(self, moment: datetime) -> bool:
        """
        Whether ``moment`` falls within this period, both bounds inclusive.

        Naive datetimes on either side are read as UTC, matching the rest of the pipeline.

        :param moment: point in time to test.
        :returns: True if the period is in force at ``moment``.
        """
        moment = as_aware_datetime(moment)

        if moment < as_aware_datetime(self.start):
            return False

        if self.end is None:
            # a period without an end never expires
            return True

        return moment <= as_aware_datetime(self.end)


class Policy(StrictIgnoringBaseModel):
    uri: str


class Coding(FhirElement):
    system: str
    version: str | None = None
    code: str
    display: str | None = None
    user_selected: bool | None = None


class CodeableConcept(FhirElement):
    coding: Annotated[list[Coding], Field(min_length=1)]
    text: str | None = None

    @property
    def codings(self) -> frozenset[tuple[str, str]]:
        """The (system, code) pairs carried by this concept."""
        return frozenset((coding.system, coding.code) for coding in self.coding)


class ConsentProvision(StrictIgnoringBaseModel):
    type: ProvisionType
    period: Period
    code: Annotated[list[CodeableConcept], Field(min_length=1)]
    #: Forbidden by the profile (0..0), and rejected rather than ignored: evaluation never reads it.
    provision: list[Any] | None = None

    @model_validator(mode="after")
    def reject_nested_provision(self):
        if self.provision:
            raise ValueError(
                "consent.provision.provision[].provision is not allowed by the MII consent profile "
                "and its permissions would not be evaluated"
            )
        return self


class RootConsentProvision(StrictIgnoringBaseModel):
    type: ProvisionType
    #: Required by every profile version; only its end became optional in profile 1.0.9.
    period: Period
    provision: list[ConsentProvision] = Field(default_factory=list)
    #: The profile constrains Consent.provision.code to 0..0; permissions belong on the sub-provisions.
    code: list[Any] | None = None

    @model_validator(mode="after")
    def reject_root_code(self):
        if self.code:
            raise ValueError(
                "consent.provision.code is not allowed by the MII consent profile; "
                "permissions belong on consent.provision.provision[].code"
            )
        return self


class Identifier(StrictIgnoringBaseModel):
    # The profile requires both on Consent.patient.identifier (1..1 each).
    system: str
    value: str


class Patient(StrictIgnoringBaseModel):
    # The profile marks both as mustSupport yet requires neither, so the patient may be identified
    # either way. A patient stating neither (e.g. display only) carries nothing this model keeps
    # and is rejected rather than silently reduced to an empty object.
    reference: str | None = None
    identifier: Identifier | None = None

    @model_validator(mode="after")
    def require_reference_or_identifier(self):
        if self.reference is None and self.identifier is None:
            raise ValueError("consent.patient must identify the patient by a reference or an identifier")
        return self


class Verification(StrictIgnoringBaseModel):
    verified: bool
    verification_date: datetime | None = None


EXPECTED_SCOPE_CODING_SYSTEM = "http://terminology.hl7.org/CodeSystem/consentscope"
EXPECTED_SCOPE_CODING_CODE = "research"
MII_BROAD_CONSENT_OID = "2.16.840.1.113883.3.1937.777.24.2.184"
#: OIDs of the broad consent versions and additional modules, introduced in package 2026.0.0.
MII_CONSENT_VERSION_MODULES_SYSTEM = (
    "https://www.medizininformatik-initiative.de/fhir/modul-consent/CodeSystem/mii-cs-consent-version-modules"
)
#: The category CodeSystem was renamed in package version 2026.0.0, but that package's own
#: examples still use the old spelling, so both are in the wild.
MII_CONSENT_CATEGORY_SYSTEMS = (
    "https://www.medizininformatik-initiative.de/fhir/modul-consent/CodeSystem/mii-cs-consent-consent_category",
    MII_CONSENT_VERSION_MODULES_SYSTEM,
)
EXPECTED_CATEGORIES: dict[str, frozenset[tuple[str, str]]] = {
    "loinc": frozenset({("http://loinc.org", "57016-8")}),
    "mii": frozenset({(s, MII_BROAD_CONSENT_OID) for s in MII_CONSENT_CATEGORY_SYSTEMS}),
}


class BroadConsentVersion(StrEnum):
    """Version of the MII broad consent document a donor was presented with."""

    V1_6D = "1.6d"
    V1_6F = "1.6f"
    V1_7_2 = "1.7.2"


class ConsentDocumentKind(StrEnum):
    """What a broad consent document declares."""

    #: Consent to the broad consent itself.
    CONSENT = "consent"
    #: Refusal of the broad consent (Ablehnung).
    REJECTION = "rejection"
    #: Withdrawal of the entire broad consent (Komplettwiderruf).
    COMPLETE_WITHDRAWAL = "complete withdrawal"
    #: Withdrawal of parts of the broad consent (Teilwiderruf).
    PARTIAL_WITHDRAWAL = "partial withdrawal"
    #: An additional module on top of the broad consent (Zusatzmodul).
    ADDITIONAL_MODULE = "additional module"


#: What each code of the version and modules CodeSystem declares. These OIDs identify the signed
#: document and appear as Consent.policy[].uri, unlike schemaVersion, which names the KDS package.
#: Unknown OIDs are reported, never rejected: a future broad consent version must not break
#: submissions.
BROAD_CONSENT_DOCUMENT_OIDS: dict[str, tuple[ConsentDocumentKind, BroadConsentVersion | None]] = {
    MII_BROAD_CONSENT_OID: (ConsentDocumentKind.CONSENT, None),
    "2.16.840.1.113883.3.1937.777.24.2.1790": (ConsentDocumentKind.CONSENT, BroadConsentVersion.V1_6D),
    "2.16.840.1.113883.3.1937.777.24.2.4053": (ConsentDocumentKind.REJECTION, BroadConsentVersion.V1_6D),
    "2.16.840.1.113883.3.1937.777.24.2.2718": (ConsentDocumentKind.COMPLETE_WITHDRAWAL, BroadConsentVersion.V1_6D),
    "2.16.840.1.113883.3.1937.777.24.2.2719": (ConsentDocumentKind.PARTIAL_WITHDRAWAL, BroadConsentVersion.V1_6D),
    "2.16.840.1.113883.3.1937.777.24.2.1791": (ConsentDocumentKind.CONSENT, BroadConsentVersion.V1_6F),
    "2.16.840.1.113883.3.1937.777.24.2.2720": (ConsentDocumentKind.COMPLETE_WITHDRAWAL, BroadConsentVersion.V1_6F),
    "2.16.840.1.113883.3.1937.777.24.2.2721": (ConsentDocumentKind.PARTIAL_WITHDRAWAL, BroadConsentVersion.V1_6F),
    "2.16.840.1.113883.3.1937.777.24.2.2079": (ConsentDocumentKind.CONSENT, BroadConsentVersion.V1_7_2),
    "2.16.840.1.113883.3.1937.777.24.2.4054": (ConsentDocumentKind.REJECTION, BroadConsentVersion.V1_7_2),
    "2.16.840.1.113883.3.1937.777.24.2.2722": (ConsentDocumentKind.COMPLETE_WITHDRAWAL, BroadConsentVersion.V1_7_2),
    "2.16.840.1.113883.3.1937.777.24.2.2723": (ConsentDocumentKind.PARTIAL_WITHDRAWAL, BroadConsentVersion.V1_7_2),
    # consent of legal guardians and of minors, all belonging to broad consent 1.7.2
    "2.16.840.1.113883.3.1937.777.24.2.3542": (ConsentDocumentKind.CONSENT, BroadConsentVersion.V1_7_2),
    "2.16.840.1.113883.3.1937.777.24.2.3543": (ConsentDocumentKind.CONSENT, BroadConsentVersion.V1_7_2),
    "2.16.840.1.113883.3.1937.777.24.2.3544": (ConsentDocumentKind.CONSENT, BroadConsentVersion.V1_7_2),
    # additional modules, not tied to a single broad consent version
    "2.16.840.1.113883.3.1937.777.24.2.4052": (ConsentDocumentKind.ADDITIONAL_MODULE, None),
    "2.16.840.1.113883.3.1937.777.24.2.4031": (ConsentDocumentKind.ADDITIONAL_MODULE, None),
    "2.16.840.1.113883.3.1937.777.24.2.4036": (ConsentDocumentKind.ADDITIONAL_MODULE, None),
    "2.16.840.1.113883.3.1937.777.24.2.4037": (ConsentDocumentKind.ADDITIONAL_MODULE, None),
    "2.16.840.1.113883.3.1937.777.24.2.4048": (ConsentDocumentKind.ADDITIONAL_MODULE, None),
}

_OID_URI_PREFIX = "urn:oid:"


def _strip_oid_prefix(uri: str) -> str:
    """Reduce a policy URI to a bare OID; the MII ships examples both with and without the prefix.

    :param uri: policy URI as submitted.
    :returns: the URI without a leading 'urn:oid:'.
    """
    return uri.removeprefix(_OID_URI_PREFIX)


class Consent(StrictIgnoringBaseModel):
    resource_type: Literal["Consent"] | None = None
    status: Status
    scope: CodeableConcept
    category: Annotated[list[CodeableConcept], Field(min_length=2)]
    patient: Patient
    date_time: datetime
    policy: Annotated[list[Policy], Field(min_length=1)]
    verification: list[Verification] | None = None
    provision: RootConsentProvision | None = None

    @model_validator(mode="after")
    def ensure_valid_scope(self):
        if len(self.scope.coding) != 1:
            raise ValueError(f"consent.scope.coding must contain only a single element, not {len(self.scope.coding)}")

        if self.scope.coding[0].system != EXPECTED_SCOPE_CODING_SYSTEM:
            raise ValueError(
                f"Expected '{EXPECTED_SCOPE_CODING_SYSTEM}' in consent.scope.coding[0].system, got '{self.scope.coding[0].system}'"
            )

        if self.scope.coding[0].code != EXPECTED_SCOPE_CODING_CODE:
            raise ValueError(
                f"Expected '{EXPECTED_SCOPE_CODING_CODE}' in consent.scope.coding[0].code, got '{self.scope.coding[0].code}'"
            )

        return self

    @model_validator(mode="after")
    def ensure_valid_category(self):
        categories_to_find = set(EXPECTED_CATEGORIES.keys())
        for i, category in enumerate(self.category):
            for expected_category_name, accepted in EXPECTED_CATEGORIES.items():
                if not (category.codings & accepted):
                    continue
                # only the loinc and mii slices are pinned to one coding; the slicing is open
                if len(category.coding) != 1:
                    raise ValueError(
                        f"consent.category[{i}] carries the {expected_category_name} category and must "
                        f"therefore contain only a single coding, not {len(category.coding)}"
                    )
                if expected_category_name not in categories_to_find:
                    raise ValueError(f"Duplicate required category in consent.category: {category}")
                categories_to_find.remove(expected_category_name)

        if categories_to_find:
            missing = ", ".join(
                f"{name} ({' or '.join(f'{system}|{code}' for system, code in sorted(EXPECTED_CATEGORIES[name]))})"
                for name in sorted(categories_to_find)
            )
            raise ValueError(f"Missing expected categories: {missing}")

        return self

    @property
    def is_in_force(self) -> bool:
        """
        Whether the consent's status marks it as currently valid.

        Consent.status is a FHIR modifier element: every status other than 'active' marks the consent
        as not in force, so its provisions must not be read as permissions.
        """
        return self.status == Status.ACTIVE

    @property
    def document_oids(self) -> frozenset[str]:
        """
        Bare OIDs identifying the signed documents, from policy URIs and from category codings in
        the version and modules CodeSystem.
        """
        from_policies = {_strip_oid_prefix(policy.uri) for policy in self.policy}
        from_categories = {
            coding.code
            for category in self.category
            for coding in category.coding
            if coding.system == MII_CONSENT_VERSION_MODULES_SYSTEM
        }
        return frozenset(from_policies | from_categories)

    @property
    def _known_documents(self) -> tuple[tuple[ConsentDocumentKind, BroadConsentVersion | None], ...]:
        """The (kind, version) pairs of those document OIDs the version and modules CodeSystem defines."""
        return tuple(
            BROAD_CONSENT_DOCUMENT_OIDS[oid] for oid in self.document_oids if oid in BROAD_CONSENT_DOCUMENT_OIDS
        )

    @property
    def broad_consent_versions(self) -> frozenset[BroadConsentVersion]:
        """Broad consent versions this consent declares. Empty if none of its OIDs is recognized."""
        return frozenset(version for _, version in self._known_documents if version is not None)

    @property
    def document_kinds(self) -> frozenset[ConsentDocumentKind]:
        """What the recognized OIDs of this consent declare."""
        return frozenset(kind for kind, _ in self._known_documents)

    @property
    def unknown_document_oids(self) -> frozenset[str]:
        """
        OIDs that are not part of the known version and modules CodeSystem, e.g. a newer broad consent.

        Policy URIs that are no OIDs at all are reported here unchanged.
        """
        return frozenset(oid for oid in self.document_oids if oid not in BROAD_CONSENT_DOCUMENT_OIDS)
