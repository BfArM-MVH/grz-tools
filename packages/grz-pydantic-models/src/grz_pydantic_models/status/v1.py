from __future__ import annotations

from datetime import date, datetime
from enum import StrEnum
from typing import Annotated

from pydantic import BaseModel, ConfigDict, StringConstraints

from grz_pydantic_models.submission.metadata import GenomicDataCenterId, SubmitterId

# TODO: re-use SubmissionId/FailureReason definitions from shared modules once to avoid local re-definitions here.
SubmissionId = Annotated[str, StringConstraints(pattern=r"^[0-9]{9}_\d{4}-\d{2}-\d{2}_[a-f0-9]{8}$")]


class FailureReason(StrEnum):
    DUPLICATE_TANG = "duplicate_tang"
    INCOMPLETE_SUBMISSION = "incomplete_submission"
    DECRYPTION_ERROR = "decryption_error"
    NETWORK_ERROR = "network_error"
    VALIDATION_ERROR = "validation_error"
    FILE_NOT_FOUND = "file_not_found"
    ENCRYPTION_ERROR = "encryption_error"
    UPLOAD_ERROR = "upload_error"
    UNKNOWN = "unknown"


class _StatusBaseModel(BaseModel):
    model_config = ConfigDict(extra="forbid", validate_assignment=True, use_enum_values=True)


class SubmissionStatusEntry(_StatusBaseModel):
    submission_id: SubmissionId
    submission_date: date | None = None
    latest_state: str | None = None
    latest_state_at: datetime | None = None
    failure_reason: FailureReason | None = None
    basic_qc_passed: bool | None = None
    detailed_qc_passed: bool | None = None
    has_detailed_qc_report: bool = False
    reported_date: date | None = None


class SubmissionStatusReport(_StatusBaseModel):
    generated_at: datetime
    grz_id: GenomicDataCenterId
    le_id: SubmitterId
    submissions: list[SubmissionStatusEntry]
