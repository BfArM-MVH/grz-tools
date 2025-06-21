from enum import StrEnum
from typing import Annotated

from pydantic import ConfigDict, Field

from ..common import StrictBaseModel


class StrictIgnoringBaseModel(StrictBaseModel):
    model_config = ConfigDict(extra="ignore")


class ProvisionType(StrEnum):
    DENY = "deny"
    PERMIT = "permit"


class Coding(StrictBaseModel):
    system: str = ""
    version: str | None = None
    code: str
    display: str | None = None
    user_selected: bool | None = None


class CodeableConcept(StrictBaseModel):
    coding: Annotated[list[Coding], Field(min_length=1)]
    text: str | None = None


class ConsentProvision(StrictIgnoringBaseModel):
    type: ProvisionType
    code: Annotated[list[CodeableConcept], Field(min_length=1)]


class RootConsentProvision(StrictIgnoringBaseModel):
    type: ProvisionType
    provision: list[ConsentProvision]


class Consent(StrictIgnoringBaseModel):
    provision: RootConsentProvision | None = None
