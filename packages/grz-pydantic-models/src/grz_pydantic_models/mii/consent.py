from datetime import date
from enum import StrEnum
from typing import Annotated

from pydantic import ConfigDict, Field, model_validator

from ..common import StrictBaseModel


class StrictIgnoringBaseModel(StrictBaseModel):
    model_config = ConfigDict(extra="ignore")


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


class Period(StrictBaseModel):
    start: date
    end: date


class Policy(StrictIgnoringBaseModel):
    uri: str


class Coding(StrictBaseModel):
    system: str
    version: str | None = None
    code: str
    display: str | None = None
    user_selected: bool | None = None


class CodeableConcept(StrictBaseModel):
    coding: Annotated[list[Coding], Field(min_length=1)]
    text: str | None = None


class ConsentProvision(StrictIgnoringBaseModel):
    type: ProvisionType
    period: Period
    code: Annotated[list[CodeableConcept], Field(min_length=1)]


class RootConsentProvision(StrictIgnoringBaseModel):
    type: ProvisionType
    provision: list[ConsentProvision]


class Patient(StrictIgnoringBaseModel):
    reference: str | None = None


EXPECTED_SCOPE = CodeableConcept(
    coding=[Coding(system="http://terminology.hl7.org/CodeSystem/consentscope", code="research")]
)
EXPECTED_CATEGORIES = (
    CodeableConcept(coding=[Coding(system="http://loinc.org", code="57016-8")]),
    CodeableConcept(
        coding=[
            Coding(
                system="https://www.medizininformatik-initiative.de/fhir/modul-consent/CodeSystem/mii-cs-consent-consent_category",
                code="2.16.840.1.113883.3.1937.777.24.2.184",
            )
        ]
    ),
)


class Consent(StrictIgnoringBaseModel):
    status: Status
    scope: CodeableConcept
    category: Annotated[list[CodeableConcept], Field(min_length=2)]
    patient: Patient
    date_time: date
    policy: Annotated[list[Policy], Field(min_length=1)]
    provision: RootConsentProvision | None = None

    @model_validator(mode="after")
    def ensure_valid_scope(self):
        if self.scope != EXPECTED_SCOPE:
            raise ValueError(
                "scope was not the expected value.\n"
                f"Expected: {EXPECTED_SCOPE.model_dump(exclude_defaults=True)}\n"
                f"Observed: {self.scope.model_dump(exclude_defaults=True)}\n"
            )
        return self

    @model_validator(mode="after")
    def ensure_valid_category(self):
        for category in EXPECTED_CATEGORIES:
            if category not in self.category:
                raise ValueError(
                    f"Expected {category.model_dump(exclude_defaults=True)} in categories but couldn't find it"
                )
        return self
