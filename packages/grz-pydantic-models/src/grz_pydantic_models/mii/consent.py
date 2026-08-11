import re
from calendar import monthrange
from datetime import UTC, date, datetime, time
from enum import StrEnum
from typing import Annotated, Any, Literal

from pydantic import BeforeValidator, ConfigDict, Field, model_validator

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


# Value regex of the FHIR R4 dateTime primitive, which allows a year, a year and month, a full date,
# or a date with a time. A time demands seconds and a timezone; ISO forms FHIR does not spell out,
# such as a space in place of the 'T' or the basic format 20200901, are not dateTimes.
# See https://hl7.org/fhir/R4/datatypes.html#dateTime
_FHIR_DATETIME = re.compile(
    # FHIR spells year zero out of the grammar; here `_covered_days` rejects it along with every
    # other date that does not exist, so a plain four-digit year is enough.
    r"(?P<year>[0-9]{4})"
    r"(-(?P<month>0[1-9]|1[0-2])"
    r"(-(?P<day>0[1-9]|[1-2][0-9]|3[0-1])"
    r"(?P<time>T([01][0-9]|2[0-3]):[0-5][0-9]:([0-5][0-9]|60)(\.[0-9]+)?"
    r"(Z|(\+|-)((0[0-9]|1[0-3]):[0-5][0-9]|14:00)))?)?)?"
)


def _covered_days(match: re.Match[str]) -> tuple[date, date]:
    """
    First and last day the matched FHIR date covers at the precision it states.

    :param match: match of ``_FHIR_DATETIME`` without a time.
    :returns: first and last day covered.
    :raises ValueError: if the date does not exist, which the regex alone does not rule out because
        it accepts every day from 01 to 31 in any month.
    """
    year = int(match["year"])

    if match["month"] is None:
        return date(year, 1, 1), date(year, 12, 31)
    month = int(match["month"])

    if match["day"] is None:
        return date(year, month, 1), date(year, month, monthrange(year, month)[1])

    day = date(year, month, int(match["day"]))
    return day, day


def _to_datetime(value: Any, *, last_moment: bool) -> Any:
    """
    Widen a FHIR dateTime without a time to the datetime bounding the span it names.

    FHIR reads such a value as the whole span it names, be that a day, a month or a year, so a
    period start begins at that span's first moment and a period end expires at its last. Values
    carrying a time are passed through untouched for pydantic to parse.

    A string that is not a FHIR dateTime is rejected rather than passed through: pydantic reads an
    all-numeric string as a Unix timestamp, which would place such a value in 1970 unremarked.

    The result is timezone-aware, reading a value naming no zone as UTC as the rest of the pipeline
    does. FHIR demands a zone alongside a time, so a naive result would serialize to a string this
    validator no longer accepts.

    :param value: raw value as submitted.
    :param last_moment: whether to widen to the last moment of the span ``value`` names rather than
        its first. Only a period end takes the last; a period start and a standalone point in time
        take the first.
    :returns: an aware datetime for values without a time, otherwise ``value`` unchanged.
    :raises ValueError: if ``value`` is neither ``None`` nor a FHIR dateTime.
    """
    time_of_day = time.max if last_moment else time.min

    if value is None:
        return value

    if isinstance(value, datetime):
        return as_aware_datetime(value)

    if isinstance(value, date):
        return datetime.combine(value, time_of_day, tzinfo=UTC)

    if not isinstance(value, str):
        raise ValueError(f"{type(value).__name__} is not a FHIR dateTime")

    match = _FHIR_DATETIME.fullmatch(value)
    if match is None:
        raise ValueError(
            f"'{value}' is not a FHIR dateTime; expected YYYY, YYYY-MM, YYYY-MM-DD or YYYY-MM-DDThh:mm:ss+zz:zz"
        )

    if match["time"] is not None:
        return value

    try:
        first, last = _covered_days(match)
    except ValueError:
        raise ValueError(f"'{value}' is not an existing date") from None

    return datetime.combine(last if last_moment else first, time_of_day, tzinfo=UTC)


def _at_first_moment(value: Any) -> Any:
    return _to_datetime(value, last_moment=False)


def _at_last_moment(value: Any) -> Any:
    return _to_datetime(value, last_moment=True)


FhirDateTime = Annotated[datetime, BeforeValidator(_at_first_moment)]
"""
A FHIR dateTime read as the first moment of the span it names.

This is the reading for every dateTime except a period end: a point in time bounds nothing, so one
stating only a year means the first moment that year can denote.
"""

FhirPeriodEnd = Annotated[datetime, BeforeValidator(_at_last_moment)]
"""
A FHIR dateTime read as the last moment of the span it names.

FHIR reads a period end without a time as the whole span it names, so an end stating only a year
stays in force until that year is over.
"""


class Period(FhirElement):
    start: FhirDateTime
    end: FhirPeriodEnd | None = None
    """
    Optional because one model serves every profile version and 1.0.9 relaxed this to 0..1. FHIR
    reads a period without an end as still running, so ResearchConsent rejects one when the declared
    schemaVersion ships a profile that still pins the end to 1..1.
    """

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
    provision: list[Any] | None = None
    """
    Forbidden by the profile (0..0), and rejected rather than ignored: evaluation never reads it.
    """

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
    period: Period
    """
    Required by every profile version; only its end became optional in profile 1.0.9.
    """

    provision: list[ConsentProvision] = Field(default_factory=list)

    code: list[Any] | None = None
    """
    The profile constrains Consent.provision.code to 0..0; permissions belong on the sub-provisions.
    """

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
    verification_date: FhirDateTime | None = None


EXPECTED_SCOPE_CODING_SYSTEM = "http://terminology.hl7.org/CodeSystem/consentscope"
EXPECTED_SCOPE_CODING_CODE = "research"
MII_BROAD_CONSENT_OID = "2.16.840.1.113883.3.1937.777.24.2.184"
# OIDs of the broad consent versions and additional modules, introduced in package 2026.0.0.
MII_CONSENT_VERSION_MODULES_SYSTEM = (
    "https://www.medizininformatik-initiative.de/fhir/modul-consent/CodeSystem/mii-cs-consent-version-modules"
)
# The category CodeSystem was renamed in package version 2026.0.0, but that package's own
# examples still use the old spelling, so both are in the wild.
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

    CONSENT = "consent"
    """Consent to the broad consent itself."""

    REJECTION = "rejection"
    """Refusal of the broad consent (Ablehnung)."""

    COMPLETE_WITHDRAWAL = "complete withdrawal"
    """Withdrawal of the entire broad consent (Komplettwiderruf)."""

    PARTIAL_WITHDRAWAL = "partial withdrawal"
    """Withdrawal of parts of the broad consent (Teilwiderruf)."""

    ADDITIONAL_MODULE = "additional module"
    """An additional module on top of the broad consent (Zusatzmodul)."""


# What each code of the version and modules CodeSystem declares. These OIDs identify the signed
# document and appear as Consent.policy[].uri, unlike schemaVersion, which names the KDS package.
# Unknown OIDs are reported, never rejected: a future broad consent version must not break
# submissions.
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
    date_time: FhirDateTime
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
