from datetime import date, datetime
from enum import StrEnum
from typing import Annotated

from pydantic import ConfigDict, Field, field_validator, model_validator

from ..common import StrictBaseModel


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


class Period(FhirElement):
    start: datetime
    #: Optional since MII consent package 2026.0.0: a period without an end never expires.
    end: datetime | None = None

    @field_validator("start", "end", mode="before")
    def date_to_datetime(cls, v):  # noqa: N805
        if isinstance(v, datetime):
            return v
        if isinstance(v, date):
            return datetime.combine(v, datetime.min.time())
        if isinstance(v, str):
            try:
                d = date.fromisoformat(v)
            except ValueError:
                return v
            return datetime.combine(d, datetime.min.time())
        return v


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


class ConsentProvision(StrictIgnoringBaseModel):
    type: ProvisionType
    period: Period
    code: Annotated[list[CodeableConcept], Field(min_length=1)]


class RootConsentProvision(StrictIgnoringBaseModel):
    type: ProvisionType
    provision: list[ConsentProvision] = Field(default_factory=list)


class Patient(StrictIgnoringBaseModel):
    reference: str | None = None


class Verification(StrictIgnoringBaseModel):
    verified: bool
    verification_date: datetime | None = None


EXPECTED_SCOPE_CODING_SYSTEM = "http://terminology.hl7.org/CodeSystem/consentscope"
EXPECTED_SCOPE_CODING_CODE = "research"
MII_BROAD_CONSENT_CATEGORY_CODE = "2.16.840.1.113883.3.1937.777.24.2.184"
#: The category CodeSystem was renamed in package version 2026.0.0, but that package's own
#: examples still use the old spelling, so both are in the wild.
MII_CONSENT_CATEGORY_SYSTEMS = (
    "https://www.medizininformatik-initiative.de/fhir/modul-consent/CodeSystem/mii-cs-consent-consent_category",
    "https://www.medizininformatik-initiative.de/fhir/modul-consent/CodeSystem/mii-cs-consent-version-modules",
)
EXPECTED_CATEGORIES: dict[str, frozenset[tuple[str, str]]] = {
    "loinc": frozenset({("http://loinc.org", "57016-8")}),
    "mii": frozenset({(s, MII_BROAD_CONSENT_CATEGORY_CODE) for s in MII_CONSENT_CATEGORY_SYSTEMS}),
}


class Consent(StrictIgnoringBaseModel):
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
            if len(category.coding) != 1:
                raise ValueError(
                    f"consent.category[{i}].coding must contain only a single element, not {len(category.coding)}"
                )
            coding = (category.coding[0].system, category.coding[0].code)
            for expected_category_name, accepted in EXPECTED_CATEGORIES.items():
                if coding in accepted:
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
