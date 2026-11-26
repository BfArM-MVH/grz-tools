"""A BfArM request to delete data, as ``grzctl delete`` reads it.

The fields follow the request form "Antrag auf Löschung von Daten in Phase 0 des genomDE
Modellvorhabens" (version 1.1). The LE fills it in and BfArM attaches it to its instruction
e-mail to the GRZ. The tanCs and the clinical data node on the form concern the KDK and are
not read.
"""

import datetime
from typing import Any

from grz_pydantic_models.submission.metadata.v1 import GenomicDataCenterId, SubmitterId, Tan
from pydantic import BaseModel, ConfigDict, Field, ValidationError, field_validator

PLACEHOLDER_MARKER = "<FILL IN"
"""Marks a value of the YAML file that is still to be filled in."""

DATE_PLACEHOLDER = "YYYY-MM-DD"
"""Stands for a date that is still to be filled in."""


class DeletionRequest(BaseModel):
    """A BfArM request to delete the genomic datasets of the given tanGs."""

    model_config = ConfigDict(extra="forbid")

    instructed_at: datetime.date
    """Date of BfArM's instruction e-mail, which is not on the form."""

    submitter_id: SubmitterId
    """Institutionskennzeichen of the LE that filed the request."""

    requester_name: str = Field(min_length=1)
    """Contact person at the LE."""

    requester_email: str = Field(pattern=r"^[^@\s]+@[^@\s]+\.[^@\s]+$")
    """E-mail address of the contact person."""

    requested_at: datetime.date
    """Date the LE sent the request."""

    reason: str = Field(min_length=1)
    """Why the data is to be deleted."""

    data_node_id: GenomicDataCenterId
    """The GRZ the request is addressed to."""

    tan_g: list[Tan] = Field(min_length=1)
    """tanGs of the genomic datasets to delete."""

    comment: str | None = None
    """Optional comment of the LE."""

    @field_validator("tan_g")
    @classmethod
    def _tan_g_unique(cls, tan_gs: list[str]) -> list[str]:
        seen: set[str] = set()
        duplicates: list[str] = []
        for tan_g in tan_gs:
            if tan_g.lower() in seen and tan_g not in duplicates:
                duplicates.append(tan_g)
            seen.add(tan_g.lower())
        if duplicates:
            raise ValueError(f"listed more than once: {', '.join(duplicates)}")
        return tan_gs


def deletion_request_problems(raw: dict[str, Any]) -> list[str]:
    """Describe every problem with a deletion request read from YAML, so all can be fixed at once.

    A value still holding a placeholder is reported as such, and not again by the model's checks.

    :param raw: The request as loaded from YAML.
    :returns: One description per problem, empty if the request is valid.
    """
    problems: list[str] = []
    unfilled: set[str] = set()
    for key, value in raw.items():
        values = value if isinstance(value, list) else [value]
        if any(isinstance(v, str) and PLACEHOLDER_MARKER in v for v in values):
            problems.append(f"{key}: still has the '{PLACEHOLDER_MARKER} ...>' placeholder — replace it")
            unfilled.add(key)
        elif isinstance(value, str) and value.strip().upper() == DATE_PLACEHOLDER:
            problems.append(f"{key}: still the '{DATE_PLACEHOLDER}' placeholder — use a real date, e.g. 2026-11-25")
            unfilled.add(key)

    try:
        DeletionRequest.model_validate(raw)
    except ValidationError as e:
        for error in e.errors():
            if error["loc"] and error["loc"][0] in unfilled:
                continue
            location = ".".join(str(part) for part in error["loc"]) or "request"
            problems.append(f"{location}: {error['msg']}")
    return problems
