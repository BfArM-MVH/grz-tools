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
    provision: list[ConsentProvision] = Field(default_factory=list)


class Patient(StrictIgnoringBaseModel):
    reference: str | None = None


EXPECTED_SCOPE_CODING_SYSTEM = "http://terminology.hl7.org/CodeSystem/consentscope"
EXPECTED_SCOPE_CODING_CODE = "research"
EXPECTED_CATEGORIES = {
    "loinc": CodeableConcept(coding=[Coding(system="http://loinc.org", code="57016-8")]),
    "mii": CodeableConcept(
        coding=[
            Coding(
                system="https://www.medizininformatik-initiative.de/fhir/modul-consent/CodeSystem/mii-cs-consent-consent_category",
                code="2.16.840.1.113883.3.1937.777.24.2.184",
            )
        ]
    ),
}


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
            category_minimal = CodeableConcept(
                coding=[Coding(system=category.coding[0].system, code=category.coding[0].code)]
            )
            for expected_category_name, expected_category in EXPECTED_CATEGORIES.items():
                if expected_category == category_minimal:
                    if expected_category_name not in categories_to_find:
                        raise ValueError(f"Duplicate required category in consent.category: {category}")
                    categories_to_find.remove(expected_category_name)

        if categories_to_find:
            raise ValueError(f"Missing expected categories: {[EXPECTED_CATEGORIES[c] for c in categories_to_find]}")

        return self
