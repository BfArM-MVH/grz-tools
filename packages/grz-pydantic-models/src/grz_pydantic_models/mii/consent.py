import re
from calendar import monthrange
from dataclasses import dataclass
from datetime import UTC, date, datetime, time
from enum import StrEnum
from typing import Annotated, Any, ClassVar, Literal

from pydantic import ConfigDict, Field, model_validator
from pydantic_core import core_schema

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


@dataclass(frozen=True)
class FhirDateTime:
    """
    A FHIR dateTime, parsed far enough to say what it names and whether FHIR permits it.

    The value keeps the spelling it was submitted with. It resolves to a moment only when asked,
    because the same value means different moments depending on where it sits: a period end covers
    the span it names to the last moment, everything else takes the first.

    The pattern matches a superset of the FHIR grammar: it admits a time that names no timezone,
    which FHIR forbids. Parsing one anyway is what keeps a document written before that rule was
    enforced readable, and :attr:`is_valid_fhir` is what tells a conformant value from such a one.
    """

    raw: str
    """The value exactly as submitted."""

    year: int
    month: int | None
    day: int | None
    states_a_time: bool
    """Whether the value names a time of day rather than only a year, a month or a day."""

    timezone: str | None
    """The stated zone, ``None`` when the value names none."""

    # A time demands seconds; ISO forms FHIR does not spell out, such as a space in place of the 'T'
    # or the basic format 20200901, are not dateTimes.
    # See https://hl7.org/fhir/R4/datatypes.html#dateTime
    _PATTERN: ClassVar[re.Pattern[str]] = re.compile(
        # FHIR spells year zero out of the grammar; here `covered_days` rejects it along with every
        # other date that does not exist, so a plain four-digit year is enough.
        r"(?P<year>[0-9]{4})"
        r"(-(?P<month>0[1-9]|1[0-2])"
        r"(-(?P<day>0[1-9]|[1-2][0-9]|3[0-1])"
        r"(?P<time>T([01][0-9]|2[0-3]):[0-5][0-9]:([0-5][0-9]|60)(\.[0-9]+)?"
        r"(?P<zone>Z|(\+|-)((0[0-9]|1[0-3]):[0-5][0-9]|14:00))?)?)?)?"
    )

    @classmethod
    def _try_parse(cls, value: Any) -> "FhirDateTime | None":
        """
        Parse ``value``, or return ``None`` when it is not shaped like a dateTime at all.

        :param value: candidate value; anything that is not a string is not a dateTime.
        :returns: the parsed value, or ``None``.
        """
        if not isinstance(value, str):
            return None

        match = cls._PATTERN.fullmatch(value)
        if match is None:
            return None

        return cls(
            raw=value,
            year=int(match["year"]),
            month=int(match["month"]) if match["month"] else None,
            day=int(match["day"]) if match["day"] else None,
            states_a_time=match["time"] is not None,
            timezone=match["zone"],
        )

    @classmethod
    def parse(cls, value: Any) -> "FhirDateTime":
        """
        Parse a value as submitted, or a date or datetime built in code.

        A date or datetime is spelled the way a document would state it, so a model built in code
        holds the same value as one parsed from JSON. A naive datetime is refused rather than read
        as UTC: guessing a zone is what this type exists to avoid, and code building a consent is
        the one place where the caller still knows which zone was meant. A plain date needs none,
        since FHIR asks for one only alongside a time.

        :param value: candidate value.
        :returns: the parsed value, which may still not be :attr:`is_valid_fhir`.
        :raises ValueError: if ``value`` is a naive datetime, or is not shaped like a dateTime.
        """
        if isinstance(value, datetime):
            if value.tzinfo is None:
                raise ValueError(
                    f"'{value.isoformat()}' names no timezone; FHIR requires one alongside a time, "
                    "so attach it, e.g. datetime(..., tzinfo=UTC)"
                )
            value = value.isoformat().replace("+00:00", "Z")
        elif isinstance(value, date):
            value = value.isoformat()

        if not isinstance(value, str):
            raise ValueError(f"{type(value).__name__} is not a FHIR dateTime")

        parsed = cls._try_parse(value)
        if parsed is None:
            raise ValueError(
                f"'{value}' is not a FHIR dateTime; "
                "expected YYYY, YYYY-MM, YYYY-MM-DD or YYYY-MM-DDThh:mm:ss with an optional +zz:zz"
            )

        # the pattern admits every day from 01 to 31 in any month, so check the date exists while
        # there is still a document to blame rather than when something later asks for a moment
        parsed.covered_days()
        return parsed

    @property
    def has_timezone(self) -> bool:
        """Whether the value names a timezone."""
        return self.timezone is not None

    @property
    def is_valid_fhir(self) -> bool:
        """
        Whether FHIR permits this value.

        A value that states a time must name its zone. One that states no time needs none, so a
        plain date is valid as it stands.
        """
        return self.has_timezone or not self.states_a_time

    def covered_days(self) -> tuple[date, date]:
        """
        First and last day this value covers at the precision it states.

        :returns: first and last day covered.
        :raises ValueError: if the date does not exist, which the pattern alone does not rule out
            because it accepts every day from 01 to 31 in any month.
        """
        try:
            if self.month is None:
                return date(self.year, 1, 1), date(self.year, 12, 31)

            if self.day is None:
                last = monthrange(self.year, self.month)[1]
                return date(self.year, self.month, 1), date(self.year, self.month, last)

            day = date(self.year, self.month, self.day)
            return day, day
        except ValueError:
            raise ValueError(f"'{self.raw}' is not an existing date") from None

    @property
    def first_moment(self) -> datetime:
        """
        The first moment this value covers, in UTC.

        A value naming a whole year, month or day begins at that span's first moment. One stating a
        time already names a moment, and is read as UTC when it names no zone, as
        ``as_aware_datetime`` does elsewhere in the pipeline.

        :raises ValueError: if the value names no existing date, or a time pydantic cannot represent.
        """
        return self._stated_moment() or datetime.combine(self.covered_days()[0], time.min, tzinfo=UTC)

    @property
    def last_moment(self) -> datetime:
        """
        The last moment this value covers, in UTC.

        A value naming a whole year, month or day expires at that span's last moment. One stating a
        time names a single moment, so it is both its first and its last.

        :raises ValueError: if the value names no existing date, or a time pydantic cannot represent.
        """
        return self._stated_moment() or datetime.combine(self.covered_days()[1], time.max, tzinfo=UTC)

    def _stated_moment(self) -> datetime | None:
        """The moment this value states outright, or ``None`` when it names a span instead."""
        if not self.states_a_time:
            return None

        try:
            moment = datetime.fromisoformat(self.raw)
        except ValueError:
            # the pattern admits a leap second, which datetime cannot represent
            raise ValueError(f"'{self.raw}' names a time that cannot be represented") from None

        return moment if moment.tzinfo is not None else moment.replace(tzinfo=UTC)

    def __str__(self) -> str:
        """The value as it would be written back into a document."""
        return self.raw

    @classmethod
    def __get_pydantic_core_schema__(cls, source_type: Any, handler: Any) -> core_schema.CoreSchema:
        """
        Validate through :meth:`parse` and serialize back to the submitted value, unchanged.

        Keeping the spelling is the point: a document round-trips exactly as written, so storing what
        was parsed cannot quietly turn a date into a midnight timestamp.

        `parse` is the whole validator rather than a step after a string schema, so every reason a
        value is refused is stated in one place and can say what was expected.
        """
        return core_schema.no_info_plain_validator_function(
            cls.parse,
            serialization=core_schema.plain_serializer_function_ser_schema(
                # always, not only in JSON mode: a python-mode dump would otherwise hand back this
                # type's own fields rather than the value, and warn while doing it
                str,
                return_schema=core_schema.str_schema(),
                when_used="always",
            ),
        )

    @classmethod
    def __get_pydantic_json_schema__(cls, schema: Any, handler: Any) -> dict[str, Any]:
        """
        Describe the value as the string it is, constrained to the FHIR grammar.

        Built rather than handed on, because a plain validator has no schema of its own to extend.
        ``format: date-time`` would be wrong: it promises a full timestamp, while FHIR also permits
        a year, a year and month, or a date.
        """
        # JSON Schema patterns are ECMA-262, which has no `(?P<name>...)`, and are unanchored, so a
        # substring would otherwise satisfy them.
        ecma = re.sub(r"\(\?P<[a-z]+>", "(?:", cls._PATTERN.pattern)
        return {"type": "string", "pattern": "^" + ecma + "$"}


class Period(FhirElement):
    start: FhirDateTime
    end: FhirDateTime | None = None
    """
    Optional because one model serves every profile version and 1.0.9 relaxed this to 0..1. FHIR
    reads a period without an end as still running, so ResearchConsent rejects one when the declared
    schemaVersion ships a profile that still pins the end to 1..1.
    """

    def contains(self, moment: datetime) -> bool:
        """
        Whether ``moment`` falls within this period, both bounds inclusive.

        FHIR reads a bound stating no time as the whole span it names, so a start begins at that
        span's first moment and an end expires at its last. Resolving the bounds here rather than
        when they were parsed keeps the submitted value intact and puts the asymmetry in view.

        Naive datetimes on either side are read as UTC, matching the rest of the pipeline.

        :param moment: point in time to test.
        :returns: True if the period is in force at ``moment``.
        """
        moment = as_aware_datetime(moment)

        if moment < self.start.first_moment:
            return False

        if self.end is None:
            # a period without an end never expires
            return True

        return moment <= self.end.last_moment


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

    def datetimes_fhir_does_not_permit(self) -> list[FhirDateTime]:
        """
        Every dateTime in this consent that states a time but names no timezone.

        The models keep a value as submitted, so a document written before the timezone was required
        still parses and can still be read. Refusing one is a submission-time rule, and this is what
        an ingest path asks rather than re-reading the raw document.

        :returns: the offending values, in the order encountered.
        """
        periods = []
        if self.provision is not None:
            periods.append(self.provision.period)
            periods.extend(provision.period for provision in self.provision.provision)

        candidates = [self.date_time]
        for period in periods:
            candidates.extend(bound for bound in (period.start, period.end) if bound is not None)
        candidates.extend(
            verification.verification_date
            for verification in self.verification or []
            if verification.verification_date is not None
        )

        return [candidate for candidate in candidates if not candidate.is_valid_fhir]

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
