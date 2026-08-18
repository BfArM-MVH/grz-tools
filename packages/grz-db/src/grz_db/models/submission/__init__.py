import base64
import calendar
import datetime
import logging
import math
import random
import re
from collections.abc import Generator, Sequence
from contextlib import contextmanager
from enum import Enum as PyEnum
from operator import attrgetter
from typing import Any, ClassVar, Literal, Optional, Protocol, Self, cast

import sqlalchemy as sa
import sqlalchemy.dialects.postgresql as sa_psql
from alembic import command as alembic_command
from alembic.config import Config as AlembicConfig
from alembic.runtime.migration import MigrationContext
from alembic.script import ScriptDirectory as AlembicScriptDirectory
from grz_pydantic_models.dates import date_to_quarter_year, quarter_date_bounds
from grz_pydantic_models.submission.metadata import (
    REDACTED_TAN,
    CoverageType,
    DiseaseType,
    GenomicDataCenterId,
    GenomicStudySubtype,
    GenomicStudyType,
    GrzSubmissionMetadata,
    LibraryType,
    Relation,
    ResearchConsentNoScopeJustification,
    SequenceSubtype,
    SequenceType,
    SubmissionType,
    SubmitterId,
    Tan,
    is_redacted_local_case_id,
)
from grz_pydantic_models.submission.metadata.v1 import Donor as MetadataDonor
from pydantic import ConfigDict, field_serializer, field_validator, model_validator
from sqlalchemy import JSON, BigInteger, Column, Enum
from sqlalchemy import func as sqlfn
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm import selectinload
from sqlmodel import DateTime, Field, Relationship, Session, SQLModel, create_engine, select

from ...common import (
    CaseInsensitiveStrEnum,
    ListableEnum,
    serialize_datetime_to_iso_z,
)
from ...errors import (
    AmbiguousCaseError,
    CaseHasLinkedSubmissionsError,
    CaseNotFoundError,
    DuplicateCaseError,
    DuplicateInitialSubmissionError,
    DuplicatePsnError,
    DuplicateSubmissionError,
    DuplicateTanGError,
    SubmissionBasicQCNotPassedError,
    SubmissionDateIsNoneError,
    SubmissionNotFoundError,
    SubmissionTypeInvalidForCaseError,
    SubmissionTypeIsNoneError,
)
from ..author import Author
from ..base import BaseSignablePayload, VerifiableLog
from .diff import (  # noqa: F401
    CaseLinkDiff,
    Diff,
    DiffState,
    DonorDiff,
    DonorsDiffCollection,
    FieldDiff,
    SubmissionChangeSet,
    SubmissionDiffCollection,
)

logger = logging.getLogger(__name__)


class OutdatedDatabaseSchemaError(Exception):
    pass


class SubmissionStateEnum(CaseInsensitiveStrEnum, ListableEnum):  # type: ignore[misc]
    """Submission state enum."""

    UPLOADING = "Uploading"
    UPLOADED = "Uploaded"
    DOWNLOADING = "Downloading"
    DOWNLOADED = "Downloaded"
    DECRYPTING = "Decrypting"
    DECRYPTED = "Decrypted"
    VALIDATING = "Validating"
    VALIDATED = "Validated"
    ENCRYPTING = "Encrypting"
    ENCRYPTED = "Encrypted"
    ARCHIVING = "Archiving"
    ARCHIVED = "Archived"
    REPORTING = "Reporting"
    REPORTED = "Reported"
    QCING = "QCing"
    QCED = "QCed"
    CLEANING = "Cleaning"
    CLEANED = "Cleaned"
    FINISHED = "Finished"
    ERROR = "Error"


class SubmissionStateFilterModeEnum(CaseInsensitiveStrEnum, ListableEnum):  # type: ignore[misc]
    """Submission state filter mode enum."""

    LATEST = "latest"
    ANY = "any"


class SemicolonSeparatedStringSet(sa.types.TypeDecorator):
    impl = sa.types.String

    cache_ok = True

    @property
    def python_type(self):
        return set

    def process_bind_param(self, value: set[str] | None, dialect: sa.engine.Dialect) -> str | None:
        if not value:
            # empty sets are stored as null to distinguish from a set of a single empty string
            return None

        for s in value:
            if ";" in s:
                raise ValueError(
                    f"Cannot safely serialize string '{s}' in a semicolon-separated set since it contains a semicolon."
                )

        # sort the set for consistent serialization behavior / deterministic output
        return ";".join(sorted(value))

    def process_result_value(self, value: str | None, dialect: sa.engine.Dialect) -> set[str] | None:
        return None if value is None else set(value.split(";"))


class SubmissionBase(SQLModel):
    """Submission base model."""

    model_config = ConfigDict(validate_assignment=True)  # type: ignore
    immutable_fields: ClassVar[set[str]] = {"id"}

    id: str
    tan_g: Tan | None = Field(default=None, unique=True, index=True, alias="tanG")
    pseudonym: str | None = Field(default=None, index=True)

    # fields from Prüfbericht
    submission_uploaded_date: datetime.date | None = None
    submission_type: SubmissionType | None = None
    submitter_id: SubmitterId | None = None
    data_node_id: GenomicDataCenterId | None = None
    coverage_type: CoverageType | None = None
    disease_type: DiseaseType | None = None
    basic_qc_passed: bool | None = None

    # fields also for Tätigkeitsbericht
    consented: bool | None = None
    detailed_qc_passed: bool | None = None
    genomic_study_type: GenomicStudyType | None = None
    genomic_study_subtype: GenomicStudySubtype | None = None

    # Database column indicating whether a submission is selected for in-depth QC (True/False) or not yet decided (None).
    selected_for_qc: bool | None = None

    # extra fields
    submission_size: int | None = Field(default=None, sa_type=BigInteger)
    submission_metadata: dict[str, Any] | None = Field(
        default=None,
        sa_column=Column(JSON().with_variant(sa_psql.JSONB, "postgresql")),
    )


# The submitter's local case ID is stored in the pseudonym column; the other redacted
# metadata field is named the same on both sides.
_METADATA_FIELD_TO_COLUMN = {"local_case_id": "pseudonym"}


class Submission(SubmissionBase, table=True):
    """Submission table model."""

    __tablename__ = "submissions"
    __table_args__ = {"extend_existing": True}

    id: str = Field(primary_key=True, index=True)

    @field_validator("id")
    @classmethod
    def validate_id_pattern(cls, v: str) -> str:
        pattern = r"^[0-9]{9}_\d{4}-\d{2}-\d{2}_[a-f0-9]{8}$"
        if not re.match(pattern, v):
            raise ValueError(f"Submission ID '{v}' does not match the required pattern.")
        return v

    # additionally constrained by the partial unique index ux_submissions_one_initial_per_case
    # (created in the cases migration): at most one QC-passed 'initial' submission per case
    case_id: int | None = Field(default=None, foreign_key="cases.id", index=True)

    states: list["SubmissionStateLog"] = Relationship(back_populates="submission")

    changes: list["ChangeRequestLog"] = Relationship(back_populates="submission")

    case: Optional["Case"] = Relationship(back_populates="submissions")

    def diff(
        self,
        other: Self,
        ignore_fields: set[str] | None = None,
    ) -> SubmissionDiffCollection:
        """Compare this submission against *other* and return all detected differences.

        :param other: The new submission state to compare against.
        :param ignore_fields: Field names to skip entirely during the comparison.
        :returns: A :class:`SubmissionDiffCollection` summarising all detected differences.
        :raises ValueError: If an immutable field has changed.
        """
        result = SubmissionDiffCollection()
        for key in other.model_fields_set - (ignore_fields or set()):
            old_value = getattr(self, key)
            new_value = getattr(other, key)

            # Ensure fields are cast to date
            if key == "submission_uploaded_date":
                if isinstance(old_value, datetime.datetime):
                    old_value = old_value.date()
                if isinstance(new_value, datetime.datetime):
                    new_value = new_value.date()

            field_diff = FieldDiff.classify_field(key, old_value, new_value)
            if key in other.immutable_fields and field_diff.diff.state != DiffState.UNCHANGED:
                raise ValueError(f"Column '{key}' is read-only and cannot be modified.")
            result.append(field_diff)
        return result

    def get_latest_state(self, filter_to_type: SubmissionStateEnum | None = None) -> Optional["SubmissionStateLog"]:
        states = filter(lambda state: state.state == filter_to_type, self.states) if filter_to_type else self.states
        states = sorted(states, key=attrgetter("timestamp"))
        return states[-1] if states else None

    @classmethod
    def from_metadata(
        cls,
        submission_id: str,
        metadata: GrzSubmissionMetadata,
        submission_uploaded_date: datetime.date,
    ) -> Self:
        """Construct a Submission populated with values derived from parsed metadata.

        Only the fields that can be sourced from metadata are set; system-managed
        fields (e.g. ``basic_qc_passed``, ``selected_for_qc``) are left at their
        defaults so that ``model_fields_set`` reliably indicates which fields to
        compare during a diff.
        """
        metadata_submission_date = metadata.submission.submission_date
        if isinstance(metadata_submission_date, datetime.datetime):
            metadata_submission_date = metadata_submission_date.date()

        if isinstance(submission_uploaded_date, datetime.datetime):
            submission_uploaded_date = submission_uploaded_date.date()

        return cls.model_validate(
            {
                "id": submission_id,
                "tanG": metadata.submission.tan_g,
                "submission_type": metadata.submission.submission_type,
                "submitter_id": metadata.submission.submitter_id,
                "coverage_type": metadata.submission.coverage_type,
                "disease_type": metadata.submission.disease_type,
                "genomic_study_type": metadata.submission.genomic_study_type,
                "genomic_study_subtype": metadata.submission.genomic_study_subtype,
                "pseudonym": metadata.submission.local_case_id,
                "data_node_id": metadata.submission.genomic_data_center_id,
                "consented": metadata.consents_to_research(date=metadata_submission_date),
                "submission_size": metadata.get_submission_size(),
                "submission_uploaded_date": submission_uploaded_date,
                "submission_metadata": metadata.to_redacted_dict(),
            }
        )

    def restore_redacted_fields(self, metadata: GrzSubmissionMetadata) -> frozenset[str]:
        """Restore *metadata*'s redacted fields from this row, in place.

        The counterpart to :meth:`from_metadata`: that builds a row from metadata, this
        fills metadata back in from a row. A metadata.json read back from an archive bucket
        carries redaction placeholders, while this row recorded the submitter's values when
        the submission was populated from the local, unredacted copy.

        :param metadata: Parsed metadata, mutated in place.
        :returns: Column names still redacted because this row holds no value for them.
            Pass these to :meth:`SubmissionDb.diff` as ``ignore_fields`` so that a
            placeholder is never written.
        """
        unrestored = metadata.restore_redacted_fields(tan_g=self.tan_g, local_case_id=self.pseudonym)
        return frozenset(_METADATA_FIELD_TO_COLUMN.get(field, field) for field in unrestored)


class Case(SQLModel, table=True):
    """Case table model.

    Groups the submissions for one patient. The authoritative identity is
    ``psn`` (the RKI pseudonym), unique once assigned. ``submitter_id`` and ``local_case_id`` are
    lookup keys used to resolve the case before a ``psn`` exists; they are deliberately *not* a
    uniqueness constraint, since a future flow may resolve a case by ``psn`` alone (with
    ``local_case_id`` absent).
    """

    # Without this a table model takes any value its annotations forbid, so a mistyped
    # `db case create` would store a submitter ID no submission can ever carry and the case
    # would sit unresolvable and unlinked. SQLModel builds rows through __setattr__, so this
    # covers construction as well as assignment. Loading is unaffected: SQLAlchemy populates
    # attributes without validating, so a row written before this stays readable.
    model_config = ConfigDict(validate_assignment=True)  # type: ignore[assignment]

    __tablename__ = "cases"
    __table_args__ = {"extend_existing": True}

    immutable_fields: ClassVar[set[str]] = {"id"}

    id: int | None = Field(default=None, primary_key=True)
    # uniqueness is enforced by the partial unique index ux_cases_psn (see the cases migration),
    # which SQLModel field arguments cannot express; declaring index=True here would advertise a
    # plain non-unique index that does not exist.
    psn: str | None = Field(default=None)
    submitter_id: SubmitterId | None = Field(default=None)
    local_case_id: str | None = Field(default=None)

    submissions: list["Submission"] = Relationship(back_populates="case")

    @classmethod
    def mutable_fields(cls) -> set[str]:
        """Field names that :meth:`SubmissionDb.modify_case` may change.

        All model fields except those in :attr:`immutable_fields`. Relationship
        attributes (e.g. ``submissions``) are not model fields, so they are
        naturally excluded.

        :returns: The set of modifiable field names.
        """
        return set(cls.model_fields.keys()) - cls.immutable_fields


class _Index(PyEnum):
    """A unique index of this schema, and how each backend says it was violated.

    A rejected write arrives as one ``IntegrityError`` whatever it broke, so it has to be
    classified: PostgreSQL names the index, SQLite names the indexed columns instead.
    """

    TAN_G = ("ix_submissions_tan_g", "submissions.tan_g")
    ONE_INITIAL = ("ux_submissions_one_initial_per_case", "submissions.case_id")
    CASE_KEY = ("ux_cases_submitter_local_case", "cases.submitter_id, cases.local_case_id")
    CASE_PSN = ("ux_cases_psn", "cases.psn")

    def __init__(self, index_name: str, sqlite_columns: str) -> None:
        self.index_name = index_name
        self.sqlite_columns = sqlite_columns


def _rejected_by(e: IntegrityError) -> "_Index | None":
    """Which of this schema's unique indices rejected the write, or ``None`` for anything else.

    PostgreSQL puts the index in the driver's diagnostics, which is exact, so that is preferred
    over reading the rendered message.
    """
    orig = getattr(e, "orig", None)
    diag = getattr(orig, "diag", None) if orig is not None else None
    named = getattr(diag, "constraint_name", None) if diag is not None else None
    if named is not None:
        return next((index for index in _Index if index.index_name == named), None)
    message = str(e)
    return next(
        (
            index
            for index in _Index
            if f"UNIQUE constraint failed: {index.sqlite_columns}" in message or index.index_name in message
        ),
        None,
    )


def _case_key_or_none(local_case_id: str | None) -> str | None:
    """Return *local_case_id* when it can key a case, else ``None``.

    Keying a case on a redaction placeholder would merge every redacted submission of a
    submitter into a single case (see :func:`is_redacted_local_case_id`).
    """
    return None if is_redacted_local_case_id(local_case_id) else local_case_id


class CaseResolver(Protocol):
    """Strategy for locating the case a submission belongs to.

    Different resolvers key on different identifiers, so the resolution rule can evolve
    (currently ``(submitter_id, local_case_id)``; later the RKI ``psn``) without a schema change:
    both keys are already columns on :class:`Case`. The seam is not the whole of such a switch,
    though. Only :meth:`SubmissionDb.assign_case` and :meth:`SubmissionDb.resolve_case` take a
    resolver; :meth:`SubmissionDb.diff` resolves through :data:`DEFAULT_CASE_RESOLVER`, so
    metadata-driven linking would want one threaded through there and a ``psn`` on
    :class:`CaseLinkDiff`.

    Implementations raise :class:`AmbiguousCaseError` rather than guess when their key matches
    more than one case.
    """

    def find_case(
        self,
        session: Session,
        *,
        submitter_id: str | None,
        local_case_id: str | None,
        psn: str | None,
    ) -> "Case | None": ...


class SubmitterLocalCaseResolver:
    """Resolve a case by ``(submitter_id, local_case_id)``. The current default.

    ``ux_cases_submitter_local_case`` makes the pair unique wherever both halves are present, so
    a second match means that index is not in place. Picking one of the matches would link a
    submission to another patient's case without saying so, and nothing in the key tells the
    matches apart, so resolution raises instead of guessing.
    """

    def find_case(
        self, session: Session, *, submitter_id: str | None, local_case_id: str | None, psn: str | None
    ) -> "Case | None":
        if not submitter_id or not local_case_id:
            return None
        matches = session.exec(
            select(Case).where(Case.submitter_id == submitter_id, Case.local_case_id == local_case_id)
        ).all()
        if len(matches) > 1:
            raise AmbiguousCaseError(submitter_id, local_case_id, [case.id for case in matches if case.id is not None])
        return matches[0] if matches else None


class PsnResolver:
    """Resolve a case by RKI pseudonym. For future PSN-based linking.

    Usable today by passing it to :meth:`SubmissionDb.assign_case` or
    :meth:`SubmissionDb.resolve_case`, but not reached from ``diff``/``commit_changes``, and
    nothing could feed it there yet: a psn comes from the tanG trade, not from the submitter's
    metadata.
    """

    def find_case(
        self, session: Session, *, submitter_id: str | None, local_case_id: str | None, psn: str | None
    ) -> "Case | None":
        if psn is None:
            return None
        return session.exec(select(Case).where(Case.psn == psn)).first()


DEFAULT_CASE_RESOLVER: CaseResolver = SubmitterLocalCaseResolver()


class FailureReasonEnum(CaseInsensitiveStrEnum, ListableEnum):  # type: ignore[misc]
    """Failure reason enum for submissions in ERROR state."""

    DUPLICATE_TANG = "duplicate_tang"
    DUPLICATE_INITIAL = "duplicate_initial"
    INCOMPLETE_SUBMISSION = "incomplete_submission"
    DECRYPTION_ERROR = "decryption_error"
    NETWORK_ERROR = "network_error"
    VALIDATION_ERROR = "validation_error"
    FILE_NOT_FOUND = "file_not_found"
    ENCRYPTION_ERROR = "encryption_error"
    UPLOAD_ERROR = "upload_error"
    UNKNOWN = "unknown"


class SubmissionStateLogBase(SQLModel):
    """
    Submission state log base model.
    Holds state information for each submission.
    Timestamped.
    Can optionally have associated JSON data.
    """

    state: SubmissionStateEnum
    data: dict[str, Any] | None = Field(default=None, sa_column=Column(JSON))
    grzctl_versions: dict[str, str] | None = Field(
        default=None,
        description="grzctl versions that created this state log (nullable for backward compatibility with old state logs)",
        sa_column=Column(JSON),
    )
    failure_reason: FailureReasonEnum | None = Field(
        default=None,
        sa_column=Column(
            Enum(FailureReasonEnum, values_callable=lambda e: [x.value for x in e]),
            nullable=True,
        ),
    )
    timestamp: datetime.datetime = Field(
        default_factory=lambda: datetime.datetime.now(datetime.UTC),
        sa_column=Column(DateTime(timezone=True), nullable=False),
    )

    model_config = ConfigDict(  # type: ignore
        populate_by_name=True,
    )

    @field_serializer("timestamp")
    def serialize_timestamp(self, ts: datetime.datetime) -> str:
        return serialize_datetime_to_iso_z(ts)


class SubmissionStateLogPayload(SubmissionStateLogBase, BaseSignablePayload):
    """
    Used to bundle data for signature calculation.
    """

    submission_id: str
    author_name: str


class SubmissionStateLog(SubmissionStateLogBase, VerifiableLog[SubmissionStateLogPayload], table=True):
    """Submission state log table model."""

    __tablename__ = "submission_states"
    __table_args__ = {"extend_existing": True}

    _payload_model_class: ClassVar = SubmissionStateLogPayload

    id: int | None = Field(default=None, primary_key=True, index=True)
    submission_id: str = Field(foreign_key="submissions.id", index=True)

    author_name: str = Field(index=True)
    signature: str

    submission: Submission | None = Relationship(back_populates="states")


class SubmissionStateLogCreate(SubmissionStateLogBase):
    """Submission state log create model."""

    submission_id: str
    author_name: str
    signature: str


class SubmissionCreate(SubmissionBase):
    """Submission create model."""

    id: str


class ChangeRequestEnum(CaseInsensitiveStrEnum, ListableEnum):  # type: ignore[misc]
    """Change request enum."""

    MODIFY = "Modify"
    DELETE = "Delete"
    TRANSFER = "Transfer"


class RequestRawContentType(CaseInsensitiveStrEnum, ListableEnum):  # type: ignore[misc]
    """Type/encoding of the raw (binary) content attached to a change request."""

    PDF = "PDF"
    PNG = "PNG"


_RAW_CONTENT_MAGIC: dict[RequestRawContentType, bytes] = {
    RequestRawContentType.PDF: b"%PDF-",
    RequestRawContentType.PNG: b"\x89PNG\r\n\x1a\n",
}


def detect_raw_content_type(content: bytes) -> RequestRawContentType | None:
    """Identify a raw attachment's type from its magic bytes.

    Content is authoritative — the file extension is intentionally ignored. Returns
    ``None`` if the bytes match no supported type. This is the same magic-byte table
    the model uses to verify ``request_raw_content`` against its declared type.
    """
    return next((content_type for content_type, magic in _RAW_CONTENT_MAGIC.items() if content.startswith(magic)), None)


_TEMPLATE_PLACEHOLDER_MARKER = "<FILL IN"


class ChangeRequestLogBase(SQLModel):
    """
    Base model for change request logs.
    Timestamped. Carries the audit trail (who requested the change and when)
    plus the request content as text and/or a binary blob. Optional ``data``
    JSON is reserved for type-specific extras.
    """

    change: ChangeRequestEnum
    # Audit fields are nullable in the schema so historical rows (predating these
    # columns) remain valid. Required-ness for new entries is enforced at the
    # application/CLI layer.
    requester_name: str | None = None
    requester_email: str | None = None
    requested_at: datetime.date | None = None
    request_email_content: str | None = Field(default=None, sa_column=Column(sa.Text))
    request_raw_content: bytes | None = Field(default=None, sa_column=Column(sa.LargeBinary))
    request_raw_content_type: RequestRawContentType | None = Field(
        default=None,
        sa_column=Column(Enum(RequestRawContentType, values_callable=lambda e: [x.name for x in e])),
    )
    data: dict[str, Any] | None = Field(default=None, sa_column=Column(JSON))
    timestamp: datetime.datetime = Field(
        default_factory=lambda: datetime.datetime.now(datetime.UTC),
        sa_column=Column(DateTime(timezone=True), nullable=False),
    )

    model_config = ConfigDict(  # type: ignore[assignment]
        populate_by_name=True,
    )

    @field_serializer("timestamp")
    def serialize_timestamp(self, ts: datetime.datetime) -> str:
        return serialize_datetime_to_iso_z(ts)

    @field_serializer("request_raw_content", when_used="json")
    def _serialize_raw_content(self, v: bytes | None) -> str | None:
        return None if v is None else base64.b64encode(v).decode("ascii")

    @field_validator("request_raw_content", mode="before")
    @classmethod
    def _decode_raw_content(cls, v: Any) -> Any:
        # Accept base64 strings (e.g. when reconstructing the payload from a model_dump
        # round-trip during signature verification) as well as raw bytes (DB / CLI path).
        if isinstance(v, str):
            return base64.b64decode(v)
        return v

    @model_validator(mode="after")
    def _reject_template_placeholders(self) -> Self:
        for name in ("requester_name", "requester_email", "request_email_content"):
            value = getattr(self, name)
            if isinstance(value, str) and _TEMPLATE_PLACEHOLDER_MARKER in value:
                raise ValueError(
                    f"Field '{name}' still contains the template placeholder "
                    f"'{_TEMPLATE_PLACEHOLDER_MARKER}…>'. Replace it with the actual value."
                )
        return self

    @model_validator(mode="after")
    def _check_raw_content_pair(self) -> Self:
        # Note: "at least one of email_content / raw_content" is enforced at the CLI layer
        # so that historical rows (with both NULL) can still be loaded from the DB.
        if (self.request_raw_content is None) != (self.request_raw_content_type is None):
            raise ValueError(
                "'request_raw_content' and 'request_raw_content_type' must be set together (or both omitted)."
            )
        if self.request_raw_content is not None and self.request_raw_content_type is not None:
            expected_magic = _RAW_CONTENT_MAGIC.get(self.request_raw_content_type)
            if expected_magic is not None and not self.request_raw_content.startswith(expected_magic):
                raise ValueError(
                    f"request_raw_content does not start with the expected "
                    f"{self.request_raw_content_type.value} magic bytes — file content does not "
                    f"match the declared type."
                )
        return self


class ChangeRequestLogPayload(ChangeRequestLogBase, BaseSignablePayload):
    """
    Used to bundle data for signature calculation.
    """

    submission_id: str
    author_name: str

    def to_bytes(self) -> bytes:
        # exclude_none keeps signatures stable across schema additions: rows signed before
        # the audit columns existed were signed without those keys, so re-verifying them
        # under the expanded model must also serialize without the (now-NULL) keys.
        payload_json = self.model_dump_json(by_alias=True, exclude_none=True)
        return payload_json.encode("utf8")


class ChangeRequestLog(ChangeRequestLogBase, VerifiableLog[ChangeRequestLogPayload], table=True):
    """Change-request log table model."""

    __tablename__ = "submission_change_requests"
    __table_args__ = (
        sa.CheckConstraint(
            "(request_raw_content IS NULL) = (request_raw_content_type IS NULL)",
            name="chk_change_request_raw_content_type_paired",
        ),
        {"extend_existing": True},
    )

    _payload_model_class: ClassVar = ChangeRequestLogPayload

    id: int | None = Field(default=None, primary_key=True, index=True)
    submission_id: str = Field(foreign_key="submissions.id", index=True)

    author_name: str = Field(index=True)
    signature: str

    submission: Submission | None = Relationship(back_populates="changes")


class ChangeRequestLogCreate(ChangeRequestLogBase):
    """Change request log create model."""

    submission_id: str
    author_name: str
    signature: str


def coerce_empty_set_to_none(value: set | None) -> set | None:
    """SemicolonSeparatedStringSet stores both empty sets and None as None."""
    return value if value else None


class Donor(SQLModel, table=True):
    """Donor database model."""

    __tablename__ = "donors"
    __table_args__ = {"extend_existing": True}

    submission_id: str = Field(foreign_key="submissions.id", primary_key=True)
    pseudonym: str = Field(primary_key=True)
    # use values_callable so enum name is explicitly stored across all
    # dialects. SQLite stores string of member name without any enforcement on
    # values because it doesn't have native Enum support, so this keeps things
    # consistent with other SQL server implementations.
    relation: Relation = Field(sa_column=Column(Enum(Relation, values_callable=lambda e: [x.name for x in e])))
    library_types: set[LibraryType] = Field(sa_column=Column(SemicolonSeparatedStringSet))
    sequence_types: set[SequenceType] = Field(sa_column=Column(SemicolonSeparatedStringSet))
    sequence_subtypes: set[SequenceSubtype] = Field(sa_column=Column(SemicolonSeparatedStringSet))
    mv_consented: bool
    research_consented: bool | None = None
    research_consent_missing_justifications: set[ResearchConsentNoScopeJustification] | None = Field(
        default=None, sa_column=Column(SemicolonSeparatedStringSet, nullable=True)
    )

    @field_validator("research_consent_missing_justifications")
    @classmethod
    def validate_and_coerce_justifications(cls, v: set | None) -> set | None:
        return coerce_empty_set_to_none(v)

    @classmethod
    def from_donor_metadata(
        cls,
        submission_id: str,
        donor: MetadataDonor,
        metadata_submission_date: datetime.date,
    ) -> Self:
        """Construct a Donor populated from parsed metadata.

        :param submission_id: Submission ID.
        :param donor: Donor metadata.
        :param metadata_submission_date: Submission date from metadata.
            Used to evaluate research consent at the time the submission was created.
        """
        return cls.model_validate(
            dict(
                submission_id=submission_id,
                pseudonym="index" if donor.relation == Relation.index_ else donor.donor_pseudonym,
                relation=Relation(donor.relation),
                library_types={datum.library_type for datum in donor.lab_data},
                sequence_types={datum.sequence_type for datum in donor.lab_data},
                sequence_subtypes={datum.sequence_subtype for datum in donor.lab_data},
                mv_consented=donor.consents_to_mv(),
                research_consented=donor.consents_to_research(date=metadata_submission_date),
                research_consent_missing_justifications={
                    consent.no_scope_justification
                    for consent in donor.research_consents
                    if consent.no_scope_justification is not None
                }
                if donor.research_consents
                else None,
            )
        )


class DetailedQCResult(SQLModel, table=True):
    """Detailed QC pipeline result model."""

    __tablename__ = "detailed_qc_results"
    __table_args__ = (
        sa.ForeignKeyConstraint(["submission_id", "pseudonym"], ["donors.submission_id", "donors.pseudonym"]),
        {"extend_existing": True},
    )

    submission_id: str = Field(primary_key=True)
    lab_datum_id: str = Field(primary_key=True)
    pseudonym: str
    timestamp: datetime.datetime = Field(
        default_factory=lambda: datetime.datetime.now(datetime.UTC),
        sa_column=Column(DateTime(timezone=True), nullable=False, primary_key=True),
    )
    sequence_type: SequenceType
    sequence_subtype: SequenceSubtype
    library_type: LibraryType
    percent_bases_above_quality_threshold_minimum_quality: float
    percent_bases_above_quality_threshold_percent: float
    percent_bases_above_quality_threshold_passed_qc: bool
    percent_bases_above_quality_threshold_percent_deviation: float
    mean_depth_of_coverage: float
    mean_depth_of_coverage_passed_qc: bool
    mean_depth_of_coverage_percent_deviation: float
    targeted_regions_min_coverage: float
    targeted_regions_above_min_coverage: float
    targeted_regions_above_min_coverage_passed_qc: bool
    targeted_regions_above_min_coverage_percent_deviation: float
    qc_workflow_version: str | None = Field(
        default=None,
        sa_column=Column(sa.String(length=64), nullable=True),
        description="QC workflow version (nullable for backward compatibility with old results)",
    )

    model_config = ConfigDict(  # type: ignore
        populate_by_name=True,
    )

    @field_serializer("timestamp")
    def serialize_timestamp(self, ts: datetime.datetime) -> str:
        return serialize_datetime_to_iso_z(ts)


class QCQueueEntry(SQLModel, table=True):
    """Queue of submissions that passed basic QC, ordered by pass timestamp."""

    __tablename__ = "qc_queue"
    __table_args__ = {"extend_existing": True}

    submission_id: str = Field(foreign_key="submissions.id", primary_key=True, index=True)
    basic_qc_passed_at: datetime.datetime = Field(
        default_factory=lambda: datetime.datetime.now(datetime.UTC),
        sa_column=Column(DateTime(timezone=True), nullable=False, index=True),
    )


class SubmissionDb:
    """
    API entrypoint for managing submissions.
    """

    def __init__(self, db_url: str, author: Author | None, debug: bool = False):
        """
        Initializes the SubmissionDb.

        Args:
            db_url: Database URL.
            debug: Whether to echo SQL statements.
        """
        self.engine = create_engine(db_url, echo=debug)
        self._author = author
        self._schema_confirmed = False

    @contextmanager
    def transaction(self, session: Session | None = None) -> Generator[Session, Any, None]:
        """Open a session that commits when its body completes, or join one already open.

        Joining is what lets a caller apply several writes as one transaction: whoever opened
        the session is the one that commits it, so a write that joins leaves that decision to
        its caller, and a rejected write discards everything done under it. A write therefore
        only ever flushes, which is also what puts it in front of the indices while the
        handler that can name what they rejected is still in scope.

        Instances stay usable after the commit, since nothing here is filled in by the
        database itself.

        :param session: A session to join, or ``None`` to open one.
        :raises OutdatedDatabaseSchemaError: if the database is not at the latest revision.
        """
        if session is not None:
            yield session
            return
        self._confirm_schema()
        with Session(self.engine, expire_on_commit=False) as own_session:
            yield own_session
            own_session.commit()

    def _get_alembic_config(self) -> AlembicConfig:
        """
        Loads the alembic configuration.

        Args:
            alembic_ini_path: Path to alembic ini file.
        """
        alembic_cfg = AlembicConfig()
        alembic_cfg.set_main_option("script_location", "grz_db:migrations")
        alembic_cfg.set_main_option("sqlalchemy.url", str(self.engine.url))
        return alembic_cfg

    def _confirm_schema(self) -> None:
        """Confirm, once, that the database is at the revision this code expects.

        Asking costs a scan of the migration directory and a connection of its own, and one
        pass over a submission opens a session per write, so the answer is kept. Nothing makes
        it stale in a way that helps: migrating is an operator action, downgrades are not
        supported.

        Only a confirmation is kept, never a failure. That is what lets ``db init`` and
        ``db upgrade`` build this object while the schema is behind, which is the one situation
        they exist for, and use it once they have fixed it.

        :raises OutdatedDatabaseSchemaError: if the database is not at the latest revision.
        """
        if self._schema_confirmed:
            return
        if not self._at_latest_schema():
            raise OutdatedDatabaseSchemaError(
                "Database not at latest schema. Please backup the database and then attempt a migration with `grzctl db upgrade`."
            )
        self._schema_confirmed = True

    def _at_latest_schema(self) -> bool:
        directory = AlembicScriptDirectory.from_config(self._get_alembic_config())
        with self.engine.connect() as connection:
            context = MigrationContext.configure(connection)
            return set(context.get_current_heads()) == set(directory.get_heads())

    def initialize_schema(self):
        """Initialize the database."""
        self.upgrade_schema()

    def upgrade_schema(self, revision: str = "head"):
        """
        Upgrades the database schema using alembic.

        Args:
            alembic_ini_path: Path to the alembic.ini file.
            revision: The Alembic revision to upgrade to (default: 'head').

        Raises:
            RuntimeError: For underlying Alembic errors.
        """
        alembic_cfg = self._get_alembic_config()
        try:
            alembic_command.upgrade(alembic_cfg, revision)
        except Exception as e:
            raise RuntimeError(f"Alembic upgrade failed: {e}") from e

    def add_submission(
        self,
        submission_id: str,
    ) -> Submission:
        """
        Adds a submission to the database.

        Args:
            submission_id: Submission ID.

        Returns:
            An instance of Submission.
        """
        with self.transaction() as session:
            existing_submission = session.get(Submission, submission_id)
            if existing_submission:
                raise DuplicateSubmissionError(submission_id)

            submission_create = SubmissionCreate(id=submission_id)
            db_submission = Submission.model_validate(submission_create)

            session.add(db_submission)
            session.flush()
            return db_submission

    @contextmanager
    def _translating_conflicts(
        self,
        session: Session,
        *,
        case_id: int | None = None,
        case_key: tuple[str | None, str | None] | None = None,
        psn: str | None = None,
    ) -> Generator[None, None, None]:
        """Turn a write this schema's indices rejected into the error that says what was rejected.

        Every value is taken before the body runs, because the rollback expires the instances
        they would otherwise be read from afterwards. That is why callers may pass an attribute
        of a pending row straight in: an argument is evaluated at the ``with``, before the flush
        that can be rejected.

        :param session: Session to roll back, and to look the blocking submission up in.
        :param case_id: Case whose one-initial slot the write targets.
        :param case_key: ``(submitter_id, local_case_id)`` the write would give a case.
        :param psn: Pseudonym the write would give a case.
        :raises DuplicateTanGError: if a submission's ``tan_g`` is already stored.
        :raises DuplicateInitialSubmissionError: if the case already has a QC-passed initial.
        :raises DuplicateCaseError: if a case already holds ``case_key``.
        :raises DuplicatePsnError: if a case already holds ``psn``.
        """
        try:
            yield
        except IntegrityError as e:
            session.rollback()
            match _rejected_by(e):
                case _Index.TAN_G:
                    raise DuplicateTanGError() from e
                case _Index.ONE_INITIAL:
                    raise self._duplicate_initial_error(session, case_id) from e
                case _Index.CASE_KEY:
                    raise DuplicateCaseError(*cast(tuple[str | None, str | None], case_key)) from e
                case _Index.CASE_PSN:
                    raise DuplicatePsnError(cast(str, psn)) from e
                case None:
                    raise
        except Exception:
            session.rollback()
            raise

    def _duplicate_initial_error(self, session: Session, case_id: int | None) -> DuplicateInitialSubmissionError:
        """Build the domain error for a write the ``ux_submissions_one_initial_per_case`` index rejected.

        The preceding rollback expires the session, so the submission that already passed basic QC
        is queried again to name it in the message.

        :param session: The session whose transaction was rolled back.
        :param case_id: The case the rejected write targeted, read before the commit. A unique
            index leaves NULLs unconstrained, so a rejected write always names a case.
        """
        qc_passed_initial = self._qc_passed_initial_of(session, case_id) if case_id is not None else None
        return DuplicateInitialSubmissionError(case_id, qc_passed_initial.id if qc_passed_initial is not None else None)

    def modify_submission(self, submission_id: str, key: str, value: Any, session: Session | None = None) -> Submission:
        """Set one column of a submission.

        :param session: Transaction to join; a fresh one is opened and committed when absent.
        """
        if key not in SubmissionBase.model_fields:
            raise ValueError(f"Unknown column key '{key}'")
        elif key in SubmissionBase.immutable_fields:
            raise ValueError(f"Column '{key}' is read-only and cannot be modified.")

        with self.transaction(session) as active_session:
            submission = active_session.get(Submission, submission_id)
            if submission is None:
                raise SubmissionNotFoundError(submission_id)

            setattr(submission, key, value)

            # The queue lookup below autoflushes the pending write, so the index can reject it
            # there rather than at the commit; both paths must reach the handler.
            with self._translating_conflicts(active_session, case_id=submission.case_id):
                if key == "basic_qc_passed":
                    # Basic QC state changed -> Align the in-depth QC queue to the new state
                    queue_entry = active_session.get(QCQueueEntry, submission_id)

                    if submission.basic_qc_passed is True and queue_entry is None:
                        # Basic QC passed -> Ensure that submission is tracked in the in-depth QC queue
                        active_session.add(QCQueueEntry(submission_id=submission_id))
                    elif submission.basic_qc_passed is not True and queue_entry is not None:
                        # Basic QC failed -> Ensure that submission is absent from the in-depth QC queue
                        active_session.delete(queue_entry)

                    # Keep selection flag aligned with failed basic QC.
                    if submission.basic_qc_passed is False:
                        submission.selected_for_qc = False
                active_session.add(submission)
                active_session.flush()
            return submission

    def update_submission(self, submission: Submission) -> Submission:
        """
        Persists changes made to a Submission object back to the database.

        :param submission: The Submission instance with updated field values.
        :return: The updated Submission instance.
        """
        with self.transaction() as session:
            db_submission = session.get(Submission, submission.id)
            if db_submission is None:
                raise SubmissionNotFoundError(submission.id)

            for field in SubmissionBase.model_fields:
                if field in SubmissionBase.immutable_fields:
                    continue
                setattr(db_submission, field, getattr(submission, field))

            session.add(db_submission)
            with self._translating_conflicts(session, case_id=db_submission.case_id):
                session.flush()
            return db_submission

    def set_selected_for_qc(self, submission_id: str, selected_for_qc: bool) -> Submission:
        value = "true" if selected_for_qc else "false"
        return self.modify_submission(submission_id, "selected_for_qc", value)

    def _submission_counts_as_selected_for_qc(self, submission: Submission) -> bool:
        if submission.selected_for_qc is True:
            return True
        return any(state.state in (SubmissionStateEnum.QCING, SubmissionStateEnum.QCED) for state in submission.states)  # type: ignore[union-attr]

    def _list_submitter_qc_candidates(
        self,
        submitter_id: SubmitterId | None,
        start_date: datetime.date,
        end_date: datetime.date,
    ) -> Sequence[Submission]:
        with self.transaction() as session:
            return session.exec(
                select(Submission)
                .options(selectinload(Submission.states))  # type: ignore[arg-type]
                .join(QCQueueEntry, QCQueueEntry.submission_id == Submission.id)  # type: ignore[arg-type]
                .where(Submission.submission_type == SubmissionType.initial)
                .where(Submission.basic_qc_passed)  # type: ignore[arg-type]
                .where(Submission.submission_uploaded_date.between(start_date, end_date))  # type: ignore[union-attr]
                .where(Submission.submitter_id == submitter_id)
                .order_by(QCQueueEntry.basic_qc_passed_at, Submission.id)  # type: ignore[arg-type]
            ).all()

    def _is_under_qc_target(
        self,
        submissions: Sequence[Submission],
        target_proportion: float,
        period_label: str,
    ) -> bool:
        total_selected = sum(map(self._submission_counts_as_selected_for_qc, submissions))
        logger.debug(
            "Total submissions selected for QC for submitter in submission's %s: %s", period_label, total_selected
        )
        if period_label == "month":
            return not total_selected

        qc_ratio = total_selected / len(submissions)
        logger.debug(f"Total submissions for submitter in submission's {period_label}: {len(submissions)}")
        logger.debug(
            f"Ratio of submissions selected for QC for submitter in submission's {period_label}: {qc_ratio:.2%}"
        )
        return qc_ratio <= target_proportion

    def _is_randomly_selected_for_qc(
        self,
        submission: Submission,
        submissions: Sequence[Submission],
        target_proportion: float,
        salt: str | None,
    ) -> bool:
        logger.debug("Randomly choosing whether to QC or not.")
        if target_proportion <= 0:
            return False

        submission_ids = [submitter_submission.id for submitter_submission in submissions]
        try:
            absolute_index = submission_ids.index(submission.id)
        except ValueError:
            # if the submission ID isn't in the quarter list, it hasn't met the requirements to be detailed QCed
            return False

        block_size = math.floor(1 / target_proportion)
        block_index = absolute_index // block_size
        submission_quarter, submission_year = date_to_quarter_year(submission.submission_uploaded_date)  # type: ignore[arg-type]
        seed = f"{submission.submitter_id}-{submission_year}-{submission_quarter}-{block_index}-{salt}"
        rng = random.Random(seed)  # noqa: S311

        target_index_in_block = rng.randint(0, block_size - 1)
        current_index_in_block = absolute_index % block_size
        return current_index_in_block == target_index_in_block

    def update_submission_state(
        self,
        submission_id: str,
        state: SubmissionStateEnum,
        data: dict | None = None,
        grzctl_versions: dict[str, str] | None = None,
        failure_reason: FailureReasonEnum | None = None,
    ) -> SubmissionStateLog:
        """
        Updates a submission's state to the specified state.

        Args:
            submission_id: Submission ID of the submission to update.
            state: New state of the submission.
            data: Optional data to attach to the update.
            grzctl_versions: Optional dictionary of grzctl dependency versions.

        Returns:
            An instance of SubmissionStateLog.
        """
        with self.transaction() as session:
            submission = session.get(Submission, submission_id)
            if not submission:
                raise SubmissionNotFoundError(submission_id)
            if not self._author:
                raise ValueError("No author defined")

            state_log_payload = SubmissionStateLogPayload(
                submission_id=submission_id,
                author_name=self._author.name,
                state=state,
                data=data,
                grzctl_versions=grzctl_versions,
                failure_reason=failure_reason,
            )
            signature = state_log_payload.sign(self._author.private_key())

            state_log_create = SubmissionStateLogCreate(**state_log_payload.model_dump(), signature=signature.hex())
            db_state_log = SubmissionStateLog.model_validate(state_log_create)
            session.add(db_state_log)
            session.flush()
            return db_state_log

    def get_donors(
        self, submission_id: str, pseudonym: str | None = None, session: Session | None = None
    ) -> tuple[Donor, ...]:
        """Retrieve all donors for a given submission, or, optionally, only for a specific pseudonym.

        :param session: Transaction to join; a fresh one is opened and committed when absent.
        """
        with self.transaction(session) as active_session:
            statement = select(Donor).where(Donor.submission_id == submission_id)
            if pseudonym is not None:
                statement = statement.where(Donor.pseudonym == pseudonym)
            donors = tuple(active_session.exec(statement).all())
        return donors

    def add_donor(self, donor: Donor, session: Session | None = None) -> Donor:
        """Add a donor to the database.

        :param session: Transaction to join; a fresh one is opened and committed when absent.
        """
        with self.transaction(session) as active_session:
            active_session.add(donor)
            active_session.flush()
            return donor

    def update_donor(self, updated_donor: Donor, session: Session | None = None) -> Donor:
        """Update a donor in the database.

        :param session: Transaction to join; a fresh one is opened and committed when absent.
        """
        with self.transaction(session) as active_session:
            statement = (
                select(Donor)
                .where(Donor.submission_id == updated_donor.submission_id)
                .where(Donor.pseudonym == updated_donor.pseudonym)
            )
            db_donor = active_session.exec(statement).first()

            if db_donor is None:
                raise RuntimeError("Cannot update a donor that doesn't yet exist in the database.")

            if db_donor == updated_donor:
                # nothing to do
                return db_donor

            for field in Donor.model_fields:
                old_value = getattr(db_donor, field)
                new_value = getattr(updated_donor, field)
                if old_value != new_value:
                    setattr(db_donor, field, new_value)

            active_session.add(db_donor)
            active_session.flush()
            return db_donor

    def delete_donor(self, donor: Donor, session: Session | None = None) -> None:
        """Delete a donor from the database.

        :param session: Transaction to join; a fresh one is opened and committed when absent.
        """
        with self.transaction(session) as active_session:
            active_session.delete(donor)
            active_session.flush()

    def get_detailed_qc_results(self, submission_id: str) -> tuple[DetailedQCResult, ...]:
        """Retrieve all detailed QC results for a given submission."""
        with self.transaction() as session:
            statement = select(DetailedQCResult).where(DetailedQCResult.submission_id == submission_id)
            results = tuple(session.exec(statement).all())
        return results

    def add_detailed_qc_result(self, result: DetailedQCResult) -> DetailedQCResult:
        """Add or update a detailed QC result to/in the database."""
        with self.transaction() as session:
            session.add(result)

            session.flush()
            return result

    def add_change_request(  # noqa: PLR0913
        self,
        submission_id: str,
        change: ChangeRequestEnum,
        *,
        requester_name: str,
        requester_email: str,
        requested_at: datetime.date,
        request_email_content: str | None = None,
        request_raw_content: bytes | None = None,
        request_raw_content_type: RequestRawContentType | None = None,
        data: dict | None = None,
    ) -> ChangeRequestLog:
        """
        Register a change request for a submission.

        Args:
            submission_id: Submission ID of the submission to register a change request for.
            change: Requested change.
            requester_name: Full name of the requester.
            requester_email: Email address of the requester.
            requested_at: Date the change was requested.
            request_email_content: Verbatim text content of the request (optional if raw content given).
            request_raw_content: Optional binary blob (e.g. PDF bytes).
            request_raw_content_type: Type of the binary blob; required iff request_raw_content is set.
            data: Optional type-specific extras.

        Returns:
            An instance of ChangeRequestLog.
        """
        with self.transaction() as session:
            submission = session.get(Submission, submission_id)
            if not submission:
                raise SubmissionNotFoundError(submission_id)
            if not self._author:
                raise ValueError("No author defined")

            change_request_log_payload = ChangeRequestLogPayload(
                submission_id=submission_id,
                author_name=self._author.name,
                change=change,
                requester_name=requester_name,
                requester_email=requester_email,
                requested_at=requested_at,
                request_email_content=request_email_content,
                request_raw_content=request_raw_content,
                request_raw_content_type=request_raw_content_type,
                data=data,
            )
            signature = change_request_log_payload.sign(self._author.private_key())

            change_request_log_create = ChangeRequestLogCreate(
                **change_request_log_payload.model_dump(), signature=signature.hex()
            )
            db_change_request_log = ChangeRequestLog.model_validate(change_request_log_create)
            session.add(db_change_request_log)
            session.flush()
            return db_change_request_log

    def get_submission(self, submission_id: str) -> Submission | None:
        """
        Retrieves a submission and its state history.

        Args:
            submission_id: Submission ID of the submission to retrieve.

        Returns:
            An instance of Submission or None.
        """
        with self.transaction() as session:
            statement = (
                select(Submission).where(Submission.id == submission_id).options(selectinload(Submission.states))  # type: ignore[arg-type]
            )
            submission = session.exec(statement).first()
            return submission

    def get_submissions(self, submission_ids: Sequence[str]) -> list[Submission | None]:
        """Fetch a specific set of submissions by their IDs in a single query.

        :param submission_ids: IDs to fetch.
        :returns: A list of the same length as *submission_ids*, where each element is the
                  matching :class:`Submission` or ``None`` when the ID does not exist.
        """
        if not submission_ids:
            return []
        with self.transaction() as session:
            statement = (
                select(Submission)
                .where(Submission.id.in_(submission_ids))  # type: ignore[attr-defined]
                .options(selectinload(Submission.states))  # type: ignore[arg-type]
            )
            found = {s.id: s for s in session.exec(statement).all()}
        return [found.get(sid) for sid in submission_ids]

    def get_case(self, case_id: int) -> Case | None:
        """Retrieve a case by its primary key.

        :param case_id: Primary key of the case.
        :returns: The :class:`Case`, or ``None`` if no case has that ID.
        """
        with self.transaction() as session:
            return session.get(Case, case_id)

    def list_cases(self) -> list[tuple[Case, int]]:
        """List all cases together with their linked-submission count.

        Cases with no linked submissions are included with a count of ``0``.

        :returns: A list of ``(case, linked_submission_count)`` tuples, ordered
            by ascending case ID.
        """
        with self.transaction() as session:
            rows = session.exec(
                select(Case, sqlfn.count(Submission.id))  # type: ignore[arg-type]
                .outerjoin(Submission, Submission.case_id == Case.id)  # type: ignore[arg-type]
                .group_by(Case.id)  # type: ignore[arg-type]
                .order_by(Case.id)  # type: ignore[arg-type]
            ).all()
            return list(rows)

    def list_submissions_for_case(self, case_id: int) -> Sequence[Submission]:
        """List all submissions linked to a case.

        Each submission is returned with its ``states`` relationship eagerly
        loaded. An unknown ``case_id`` yields an empty result rather than an error.

        :param case_id: Primary key of the case whose submissions to list.
        :returns: Submissions linked to the case, ordered by ascending
            submission ID.
        """
        with self.transaction() as session:
            return session.exec(
                select(Submission)
                .where(Submission.case_id == case_id)
                .options(selectinload(Submission.states))  # type: ignore[arg-type]
                .order_by(Submission.id)
            ).all()

    @staticmethod
    def _qc_passed_initial_of(session: Session, case_id: int) -> Submission | None:
        """Return the case's QC-passed ``initial`` submission, if any."""
        return session.exec(
            select(Submission).where(
                Submission.case_id == case_id,
                Submission.submission_type == SubmissionType.initial,
                Submission.basic_qc_passed.is_(True),  # type: ignore[union-attr]
            )
        ).first()

    def get_initial_submission(self, case_id: int) -> Submission | None:
        """Return the ``initial`` submission of a case that passed basic QC.

        Several competing initial submissions may be linked while pending basic QC; whichever
        passes basic QC first becomes the case's QC-passed initial submission (enforced by a
        partial unique index), so the result is unambiguous.

        :param case_id: Primary key of the case.
        :returns: The case's QC-passed ``initial`` submission, or ``None`` if no linked initial
            submission has passed basic QC yet.
        """
        with self.transaction() as session:
            return self._qc_passed_initial_of(session, case_id)

    def create_case(
        self, submitter_id: str | None = None, local_case_id: str | None = None, psn: str | None = None
    ) -> Case:
        """Create a case from the given identifiers.

        The identifiers are the case's resolution keys: ``(submitter_id,
        local_case_id)`` locate a case before a ``psn`` is assigned, while
        ``psn`` (the RKI pseudonym) is the authoritative identity and must be
        unique when present. All three are optional.

        :param submitter_id: Submitter identifier resolution key.
        :param local_case_id: Submitter-local case identifier resolution key.
        :param psn: RKI pseudonym; must be unique across cases when set.
        :returns: The newly created :class:`Case`, with its database-assigned ``id`` populated.
        :raises DuplicateCaseError: if a case already has this ``(submitter_id, local_case_id)``.
        :raises DuplicatePsnError: if ``psn`` is already assigned to another case.
        """
        with self.transaction() as session:
            case = Case(submitter_id=submitter_id, local_case_id=local_case_id, psn=psn)
            session.add(case)
            # Both keys are enforced by an index, so asking first would only add a query that
            # a concurrent write can invalidate between the answer and the insert.
            with self._translating_conflicts(session, case_key=(submitter_id, local_case_id), psn=psn):
                session.flush()
            return case

    def modify_case(self, case_id: int, key: str, value: Any) -> Case:
        """Modify a single mutable field of a case.

        Only mutable columns may be changed; ``id`` (and any other field listed
        in :attr:`Case.immutable_fields`) is read-only.

        :param case_id: Primary key of the case to modify.
        :param key: Name of the column to set (e.g. ``"psn"``, ``"submitter_id"``,
            ``"local_case_id"``).
        :param value: New value for the column.
        :returns: The updated :class:`Case`, reloaded from the database.
        :raises ValueError: if ``key`` is not a column of ``cases`` or is read-only.
        :raises CaseNotFoundError: if no case has the given ``case_id``.
        :raises DuplicateCaseError: if the change would give two cases the same
            ``(submitter_id, local_case_id)``.
        :raises DuplicatePsnError: if ``key`` is ``"psn"`` and ``value`` is
            already assigned to another case.
        """
        if key in Case.immutable_fields:
            raise ValueError(f"Column '{key}' is read-only and cannot be modified.")
        if key not in Case.mutable_fields():
            raise ValueError(f"Unknown column key '{key}'")
        with self.transaction() as session:
            case = session.get(Case, case_id)
            if case is None:
                raise CaseNotFoundError(case_id)
            setattr(case, key, value)
            session.add(case)
            # Which index rejected the write says what happened; the key being edited only says
            # what was attempted, and a NOT NULL or check violation would wear the wrong name.
            with self._translating_conflicts(session, case_key=(case.submitter_id, case.local_case_id), psn=str(value)):
                session.flush()
            return case

    def delete_case(self, case_id: int) -> None:
        """Delete a case, refusing when submissions are still linked to it.

        :param case_id: Primary key of the case to delete.
        :raises CaseNotFoundError: if no case has the given ``case_id``.
        :raises CaseHasLinkedSubmissionsError: if one or more submissions are
            still linked to the case.
        """
        with self.transaction() as session:
            case = session.get(Case, case_id)
            if case is None:
                raise CaseNotFoundError(case_id)
            count = session.exec(
                select(sqlfn.count()).select_from(Submission).where(Submission.case_id == case_id)  # type: ignore[arg-type]
            ).one()
            if count:
                raise CaseHasLinkedSubmissionsError(case_id, count)
            session.delete(case)
            session.flush()

    def set_submission_case(self, submission_id: str, case_id: int) -> Submission:
        """Link (or relink) a submission to an existing case.

        A case may have at most one ``initial`` submission that passed basic QC. A partial unique
        index enforces this at commit time, so a link that would break it raises instead.

        :param submission_id: ID of the submission to link.
        :param case_id: Primary key of the target case.
        :returns: The updated :class:`Submission`, reloaded from the database.
        :raises SubmissionNotFoundError: if no submission has the given
            ``submission_id``.
        :raises SubmissionTypeInvalidForCaseError: if the submission is a ``test``
            submission, which is never case-tracked.
        :raises CaseNotFoundError: if no case has the given ``case_id``.
        :raises DuplicateInitialSubmissionError: if linking would give the case a second
            QC-passed ``initial`` submission.
        """
        with self.transaction() as session:
            submission = session.get(Submission, submission_id)
            if submission is None:
                raise SubmissionNotFoundError(submission_id)
            self._assert_case_trackable(submission_id, submission.submission_type)
            case = session.get(Case, case_id)
            if case is None:
                raise CaseNotFoundError(case_id)
            submission.case_id = case_id
            session.add(submission)
            with self._translating_conflicts(session, case_id=case_id):
                session.flush()
            return submission

    @staticmethod
    def _assert_case_trackable(submission_id: str, submission_type: SubmissionType | None) -> None:
        """Refuse a submission that cases do not cover, whichever door it arrived at.

        Excluding ``test`` submissions is the only rule cases place on a submission's type, so
        it has to hold for a link resolved from metadata and for one an operator names directly.

        An unknown type is refused too, since it may yet turn out to be ``test``: a row that
        has not been populated carries no type, and linking it would decide the question
        before the metadata answers it. Populate resolves the link itself, so nothing is lost.

        :raises SubmissionTypeInvalidForCaseError: if the submission is a ``test`` submission,
            or if its type is not known yet.
        """
        if submission_type is None:
            raise SubmissionTypeInvalidForCaseError(
                f"submission '{submission_id}' has no submission type yet, so it cannot be "
                "case-tracked; populate it first."
            )
        if submission_type == SubmissionType.test:
            raise SubmissionTypeInvalidForCaseError(
                f"'test' submission '{submission_id}' is not case-tracked; cases group clinical submissions only."
            )

    def _find_case_for_link(  # noqa: PLR0913
        self,
        session: Session,
        submission_id: str,
        *,
        submitter_id: str | None,
        local_case_id: str | None,
        psn: str | None,
        submission_type: SubmissionType,
        resolver: CaseResolver,
    ) -> tuple[Submission, "Case | None"]:
        """Locate the case a submission may link to and validate the link. Performs no writes.

        :returns: The submission and the existing matching :class:`Case`, or ``None`` in place
            of the case if none matches and creating a new one is allowed.
        :raises SubmissionNotFoundError: if no submission has the given
            ``submission_id``.
        :raises SubmissionTypeInvalidForCaseError: if the submission is a ``test``
            submission, which is never case-tracked.
        :raises AmbiguousCaseError: if the resolution key matches more than one case.
        """
        submission = session.get(Submission, submission_id)
        if submission is None:
            raise SubmissionNotFoundError(submission_id)

        self._assert_case_trackable(submission_id, submission_type)

        # No further rule by submission type. A case groups whatever arrived for one patient, so a
        # followup whose initial never reached this GRZ opens the case itself: the submission is
        # real, it is reported to BfArM on its own tanG, and refusing to link it would only keep it
        # out of the record. Competing initial submissions coexist here too; the partial unique
        # index is the sole case invariant, bounding how many may *pass* basic QC (see
        # ``ux_submissions_one_initial_per_case``, classified by :func:`_rejected_by`).
        return submission, resolver.find_case(session, submitter_id=submitter_id, local_case_id=local_case_id, psn=psn)

    def assign_case(  # noqa: PLR0913
        self,
        submission_id: str,
        *,
        submitter_id: str | None = None,
        local_case_id: str | None = None,
        psn: str | None = None,
        submission_type: SubmissionType,
        resolver: CaseResolver | None = None,
        session: Session | None = None,
    ) -> Case:
        """Resolve or create the case for a submission and link it.

        A case may have at most one ``initial`` submission that passed basic QC. Once one passes,
        no other ``initial`` submission of that case may pass.

        The *resolver* strategy decides how an existing case is located (default:
        ``(submitter_id, local_case_id)``). A brand-new case requires an ``initial`` submission.
        Competing initial submissions may share a case while pending basic QC, since the data
        alone cannot tell a re-upload from a duplicate; only one of them may ever pass basic QC.
        A non-initial submission requires the case to have an initial submission that has not
        failed basic QC.
        ``test`` submissions are never case-tracked and are rejected. Re-running for an
        already-linked submission makes no further change (safe to repeat).

        :param submission_id: ID of the submission to assign.
        :param submitter_id: Submitter identifier passed to the resolver and
            stored on a newly created case.
        :param local_case_id: Submitter-local case identifier passed to the
            resolver and stored on a newly created case.
        :param psn: RKI pseudonym passed to the resolver and stored on a newly
            created case.
        :param submission_type: Type of the submission. Any type may open a case,
            since a patient's first submission to reach this GRZ need not be the
            ``initial`` one; only ``test`` is rejected (never case-tracked).
        :param resolver: Strategy used to locate an existing case; defaults to
            :data:`DEFAULT_CASE_RESOLVER` (:class:`SubmitterLocalCaseResolver`).
        :returns: The resolved or newly created :class:`Case` the submission is
            now linked to.
        :raises SubmissionNotFoundError: if no submission has the given
            ``submission_id``.
        :raises SubmissionTypeInvalidForCaseError: if the submission is a ``test``
            submission, which is never case-tracked.
        :raises DuplicateInitialSubmissionError: if the link would give the case a
            second QC-passed initial submission.
        :raises AmbiguousCaseError: if the resolution key matches more than one case.
        """
        resolver = resolver or DEFAULT_CASE_RESOLVER
        with self.transaction(session) as active_session:
            submission, case = self._find_case_for_link(
                active_session,
                submission_id,
                submitter_id=submitter_id,
                local_case_id=local_case_id,
                psn=psn,
                submission_type=submission_type,
                resolver=resolver,
            )
            if case is None:
                # The insert reaches the database at the flush, so that is inside the block: it
                # is where a case another process created in the meantime rejects ours.
                try:
                    with self._translating_conflicts(active_session, case_key=(submitter_id, local_case_id), psn=psn):
                        case = Case(submitter_id=submitter_id, local_case_id=local_case_id, psn=psn)
                        active_session.add(case)
                        active_session.flush()
                except DuplicateCaseError:
                    # Another process created this case between our lookup and our insert. Its row
                    # is the case now, so join that one: losing the race is not a failure, and
                    # letting it surface would fail a submission that did nothing wrong. A caller's
                    # transaction loses nothing to the rollback, since the link below is the first
                    # thing this applies.
                    case = resolver.find_case(
                        active_session, submitter_id=submitter_id, local_case_id=local_case_id, psn=psn
                    )
                    if case is None:  # pragma: no cover - the row that rejected our insert must exist
                        raise
                    # The rollback expired the submission the lookup handed back.
                    submission = cast(Submission, active_session.get(Submission, submission_id))

            # Linking is a write in its own right: whether the case was found or created by
            # whoever won the race to create it, it may already hold the initial slot.
            submission.case_id = case.id
            active_session.add(submission)
            with self._translating_conflicts(active_session, case_id=case.id):
                active_session.flush()
            return case

    def resolve_case(  # noqa: PLR0913
        self,
        submission_id: str,
        *,
        submitter_id: str | None = None,
        local_case_id: str | None = None,
        psn: str | None = None,
        submission_type: SubmissionType,
        resolver: CaseResolver | None = None,
        session: Session | None = None,
    ) -> "Case | None":
        """Preview the case :meth:`assign_case` would link, without writing.

        Runs the same lookup and validation as :meth:`assign_case` but leaves
        the database unchanged.

        :param submission_id: ID of the submission to resolve a case for.
        :param submitter_id: Submitter identifier passed to the resolver.
        :param local_case_id: Submitter-local case identifier passed to the
            resolver.
        :param psn: RKI pseudonym passed to the resolver.
        :param submission_type: Type of the submission. Any type may open a case,
            since a patient's first submission to reach this GRZ need not be the
            ``initial`` one; only ``test`` is rejected (never case-tracked).
        :param resolver: Strategy used to locate an existing case; defaults to
            :data:`DEFAULT_CASE_RESOLVER` (:class:`SubmitterLocalCaseResolver`).
        :returns: The existing :class:`Case` the submission would be linked to,
            or ``None`` if :meth:`assign_case` would create a new case.
        :raises SubmissionNotFoundError: if no submission has the given
            ``submission_id``.
        :raises SubmissionTypeInvalidForCaseError: if the submission is a ``test``
            submission, which is never case-tracked.
        :raises AmbiguousCaseError: if the resolution key matches more than one case.
        """
        resolver = resolver or DEFAULT_CASE_RESOLVER
        with self.transaction(session) as active_session:
            _submission, case = self._find_case_for_link(
                active_session,
                submission_id,
                submitter_id=submitter_id,
                local_case_id=local_case_id,
                psn=psn,
                submission_type=submission_type,
                resolver=resolver,
            )
            return case

    def assert_no_duplicate_initial(  # noqa: PLR0913
        self,
        submission_id: str,
        *,
        submitter_id: str | None = None,
        local_case_id: str | None = None,
        psn: str | None = None,
        submission_type: SubmissionType,
        resolver: CaseResolver | None = None,
    ) -> None:
        """Ask, before writing anything, whether this would be a case's second QC-passed initial.

        A case may have at most one ``initial`` submission that passed basic QC, and
        ``ux_submissions_one_initial_per_case`` is what enforces it. Linking is deliberately
        permissive, so :meth:`resolve_case` never flags a duplicate and the rejection lands only
        when a second initial submission tries to *pass* basic QC. Finding out that late means
        having validated a submission that cannot be accepted, so a caller about to spend that
        effort can ask here instead.

        The case and its QC-passed initial submission are read in one transaction, which is what
        makes the answer one answer: asked separately, a competing initial can pass basic QC in
        between and this would report a submission as clear that the index is about to reject.

        Answering costs two queries, so callers that are going to write anyway should let the
        index speak rather than ask twice.

        :param submission_id: ID of the submission about to be validated. A case whose QC-passed
            initial submission *is* this one is not a duplicate.
        :param submitter_id: Submitter identifier passed to the resolver.
        :param local_case_id: Submitter-local case identifier passed to the resolver. A redaction
            placeholder cannot key a case (see :func:`_case_key_or_none`), so it is left alone.
        :param psn: RKI pseudonym passed to the resolver.
        :param submission_type: Type of the submission. Only ``initial`` can break this rule.
        :param resolver: Strategy used to locate the case; defaults to
            :data:`DEFAULT_CASE_RESOLVER`.
        :raises DuplicateInitialSubmissionError: if the case already has a different ``initial``
            submission that passed basic QC.
        :raises SubmissionNotFoundError: if no submission has the given ``submission_id``.
        :raises AmbiguousCaseError: if the resolution key matches more than one case.
        """
        if submission_type != SubmissionType.initial or _case_key_or_none(local_case_id) is None:
            return

        with self.transaction() as session:
            case = self.resolve_case(
                submission_id,
                submitter_id=submitter_id,
                local_case_id=local_case_id,
                psn=psn,
                submission_type=submission_type,
                resolver=resolver,
                session=session,
            )
            if case is None or case.id is None:
                return
            qc_passed_initial = self._qc_passed_initial_of(session, case.id)
            if qc_passed_initial is not None and qc_passed_initial.id != submission_id:
                raise DuplicateInitialSubmissionError(case.id, qc_passed_initial.id)

    def list_submissions(
        self,
        limit: int | None,
        state_filters: Sequence[SubmissionStateEnum] | None = None,
        state_filter_mode: SubmissionStateFilterModeEnum = SubmissionStateFilterModeEnum.LATEST,
    ) -> Sequence[Submission]:
        """
        Lists all submissions in the database.

        Returns:
            A list of all submissions in the database. Ordered by latest
            submission state timestamp if not null, otherwise use submission
            date, with submissions missing both of these sorting first.
        """
        with self.transaction() as session:
            latest_state_per_submission = (
                select(
                    SubmissionStateLog.submission_id.label("submission_id"),  # type: ignore[attr-defined]
                    sqlfn.max(SubmissionStateLog.timestamp).label("timestamp"),
                )
                .group_by(SubmissionStateLog.submission_id)
                .subquery("latest_state_per_submission")
            )
            statement = (
                select(Submission)
                .options(selectinload(Submission.states))  # type: ignore[arg-type]
                .join(
                    latest_state_per_submission,
                    Submission.id == latest_state_per_submission.c.submission_id,  # type: ignore[arg-type]
                    isouter=True,
                )
                .order_by(
                    sqlfn.coalesce(latest_state_per_submission.c.timestamp, Submission.submission_uploaded_date)
                    .desc()
                    .nulls_first()
                )
            )
            state_filter_values = tuple(state_filters or ())
            if state_filter_values:
                if state_filter_mode == SubmissionStateFilterModeEnum.LATEST:
                    latest_state = (
                        select(SubmissionStateLog.state)
                        .where(SubmissionStateLog.submission_id == Submission.id)
                        .order_by(
                            SubmissionStateLog.timestamp.desc(),  # type: ignore[attr-defined]
                            SubmissionStateLog.id.desc(),  # type: ignore[union-attr]
                        )  # tie-breaker
                        .limit(1)
                        .scalar_subquery()
                    )
                    statement = statement.where(latest_state.in_(state_filter_values))
                elif state_filter_mode == SubmissionStateFilterModeEnum.ANY:
                    statement = statement.where(
                        sa.exists(
                            select(1).where(
                                SubmissionStateLog.submission_id == Submission.id,
                                SubmissionStateLog.state.in_(state_filter_values),  # type: ignore[attr-defined]
                            )
                        )
                    )
                else:
                    raise ValueError(f"Unknown state_filter_mode '{state_filter_mode}'")

            if limit is not None:
                statement = statement.limit(limit)

            submissions = session.exec(statement).all()
            return submissions

    def list_processed_between(self, start: datetime.date, end: datetime.date) -> Sequence[Submission]:
        """
        Lists all submissions processed between the given start and end dates, inclusive.
        Processed is defined as either reported (Prüfbericht submitted) or detailed QC finished.
        """
        with self.transaction() as session:
            reported_within_window = (
                select(SubmissionStateLog.submission_id)
                .where(SubmissionStateLog.state.in_([SubmissionStateEnum.REPORTED, SubmissionStateEnum.QCED]))  # type: ignore[attr-defined]
                .where(SubmissionStateLog.timestamp.between(start, end))  # type: ignore[attr-defined]
                .distinct()
                .subquery()
            )
            statement = (
                select(Submission)
                .options(selectinload(Submission.states))  # type: ignore[arg-type]
                .join(reported_within_window, Submission.id == reported_within_window.c.submission_id)  # type: ignore[arg-type]
            )
            submissions = session.exec(statement).all()
            return submissions

    def list_change_requests(self) -> Sequence[Submission]:
        """
        Lists all submissions in the database.

        Returns:
            A list of all submissions in the database, ordered by their ID.
        """
        with self.transaction() as session:
            statement = (
                select(Submission)
                .where(Submission.changes.any())  # type: ignore[attr-defined]
                .options(selectinload(Submission.changes))  # type: ignore[arg-type]
                .order_by(Submission.id)
            )
            change_requests = session.exec(statement).all()
            return change_requests

    def should_qc(self, submission_id: str, target_percentage: float, salt: str | None) -> bool:  # noqa: C901
        """
        Determines whether or not a submission should go through detailed QC or not.
        """
        target_proportion = target_percentage / 100.0
        submission = self.get_submission(submission_id)

        if submission is None:
            raise SubmissionNotFoundError(submission_id)
        submission_date = submission.submission_uploaded_date
        if submission_date is None:
            raise SubmissionDateIsNoneError()
        submission_type = submission.submission_type
        if submission_type is None:
            raise SubmissionTypeIsNoneError()
        if submission_type != SubmissionType.initial:
            # only initial submissions matter for detailed QC selection
            return False
        if submission.basic_qc_passed is not True:
            # only submissions that passed basic QC are eligible for detailed QC
            raise SubmissionBasicQCNotPassedError(submission_id)
        if submission.selected_for_qc is True:
            return True
        if submission.selected_for_qc is False:
            return False

        submission_month = submission_date.month
        submission_quarter, submission_year = date_to_quarter_year(submission_date)
        submission_quarter_start, submission_quarter_end = quarter_date_bounds(
            quarter=submission_quarter, year=submission_year
        )
        _, days_in_submission_month = calendar.monthrange(submission_year, submission_month)
        should_select = False

        # yes if none QCed/QCing from submitter yet for the submission month
        submitter_submissions_month = self._list_submitter_qc_candidates(
            submitter_id=submission.submitter_id,
            start_date=datetime.date(year=submission_year, month=submission_month, day=1),
            end_date=datetime.date(year=submission_year, month=submission_month, day=days_in_submission_month),
        )
        if self._is_under_qc_target(submitter_submissions_month, target_proportion, period_label="month"):
            should_select = True

        # yes if we are under target percentage for submitter for the submission's quarter
        if not should_select:
            submitter_submissions_quarter = self._list_submitter_qc_candidates(
                submitter_id=submission.submitter_id,
                start_date=submission_quarter_start,
                end_date=submission_quarter_end,
            )
            if self._is_under_qc_target(
                submitter_submissions_quarter,
                target_proportion,
                period_label="quarter",
            ):
                should_select = True

        # randomly, but reproducibly, select submissions for a given submitter, quarter, block, and salt
        if not should_select:
            should_select = self._is_randomly_selected_for_qc(
                submission=submission,
                submissions=submitter_submissions_quarter,
                target_proportion=target_proportion,
                salt=salt,
            )

        self.set_selected_for_qc(submission_id, should_select)
        return should_select

    @staticmethod
    def _diff_metadata(
        current_submission: Submission,
        metadata: GrzSubmissionMetadata,
        submission_uploaded_date: datetime.date,
        ignore_fields: set[str] | None = None,
    ) -> SubmissionDiffCollection:
        """Compare a submission's current database state against fresh metadata.

        :param current_submission: The submission's row, as :meth:`diff` read it.
        :param metadata: Parsed metadata from the submission's ``metadata.json``.
        :param submission_uploaded_date: Date of submission upload finish
        :param ignore_fields: Field names to skip entirely during the comparison.
        :returns: A :class:`SubmissionDiffCollection` instance summarising all detected
            differences.
        """
        if ignore_fields is None:
            ignore_fields = set()

        if isinstance(submission_uploaded_date, datetime.datetime):
            submission_uploaded_date = submission_uploaded_date.date()

        new_submission = Submission.from_metadata(current_submission.id, metadata, submission_uploaded_date)

        return current_submission.diff(new_submission, ignore_fields)

    def _diff_case_link(
        self,
        session: Session,
        current_submission: Submission,
        metadata: GrzSubmissionMetadata,
        ignore_fields: set[str],
    ) -> CaseLinkDiff | None:
        """Resolve whether committing *metadata* would change the submission's case link.

        Case resolution is skipped when ``"case_id"`` is ignored, for ``test`` submissions
        (never case-tracked), and when ``local_case_id`` is still a redaction placeholder
        (see :func:`_case_key_or_none`); call :meth:`Submission.restore_redacted_fields` first when
        *metadata* was read back from an archive bucket.

        :param session: The transaction :meth:`diff` runs in.
        :param current_submission: The submission's row, as :meth:`diff` read it.
        :param metadata: Parsed metadata from the submission's ``metadata.json``.
        :param ignore_fields: Field names skipped during the comparison.
        :returns: A :class:`CaseLinkDiff` if the case assignment would change, else ``None``.
        :raises AmbiguousCaseError: if the metadata's case key matches more than one case.
        """
        if (
            "case_id" in ignore_fields
            or metadata.submission.submission_type == SubmissionType.test
            or _case_key_or_none(metadata.submission.local_case_id) is None
        ):
            return None

        # Resolution re-reads the submission, which the session hands back from its identity
        # map: *current_submission* is held here, so the row is not fetched a second time.
        resolved = self.resolve_case(
            current_submission.id,
            submitter_id=metadata.submission.submitter_id,
            local_case_id=metadata.submission.local_case_id,
            submission_type=metadata.submission.submission_type,
            session=session,
        )
        resolved_id = resolved.id if resolved is not None else None
        if resolved is not None and resolved_id == current_submission.case_id:
            return None
        return CaseLinkDiff(
            before=current_submission.case_id,
            after=resolved_id,
            submitter_id=metadata.submission.submitter_id,
            local_case_id=metadata.submission.local_case_id,
            submission_type=metadata.submission.submission_type,
        )

    def _diff_donors(
        self,
        session: Session,
        submission_id: str,
        metadata: GrzSubmissionMetadata,
    ) -> DonorsDiffCollection:
        """Diff all donors in *metadata* against the current database state.

        :param session: The transaction :meth:`diff` runs in.
        :param submission_id: Submission ID to look up donors for.
        :param metadata: Parsed metadata from the submission's ``metadata.json``.
        :returns: A fully populated :class:`DonorDiff`.
        """
        metadata_submission_date = metadata.submission.submission_date
        if isinstance(metadata_submission_date, datetime.datetime):
            metadata_submission_date = metadata_submission_date.date()

        donors_in_db_submission = {
            donor.pseudonym: donor for donor in self.get_donors(submission_id=submission_id, session=session)
        }
        donors_in_metadata = {
            (d := Donor.from_donor_metadata(submission_id, donor, metadata_submission_date)).pseudonym: d
            for donor in metadata.donors
        }

        result = DonorsDiffCollection()

        for pseudonym in donors_in_db_submission.keys() | donors_in_metadata.keys():
            donor_before = donors_in_db_submission.get(pseudonym)
            donor_after = donors_in_metadata.get(pseudonym)

            diff = DonorDiff.classify(donor_before, donor_after)
            result.append(diff)

        return result

    def diff(
        self,
        submission_id: str,
        metadata: GrzSubmissionMetadata,
        submission_uploaded_date: datetime.date | None,
        ignore_fields: set[str] | None = None,
    ) -> SubmissionChangeSet:
        """
        Collects everything that committing *metadata* would change for a submission.

        This function compares the provided metadata and donors data with the existing
        state in the system for a specific submission ID. The returned change set carries
        the submission-level field diffs, the donor diffs, and a pending case link (see
        :attr:`SubmissionChangeSet.case_link`) when the metadata's case key resolves to a
        different case than currently linked; :meth:`commit_changes` applies the link via
        :meth:`assign_case`. ``test`` submissions are never case-tracked.

        :param submission_id: The unique identifier of the submission to be compared.
        :param metadata: The metadata of the submission to compare against.
        :param submission_uploaded_date: The date when the submission process was finished.
            If None, the field will not be included in the comparison.
        :param ignore_fields: Optional set of field names to be ignored during the metadata
            comparison. ``"case_id"`` skips case-link resolution.
        :raises SubmissionNotFoundError: if no submission has the given ``submission_id``.
        :returns: A :class:`SubmissionChangeSet` with all detected differences. A case that
            cannot be resolved is reported in
            :attr:`SubmissionChangeSet.case_link_error` rather than raised, since it says
            nothing about the other diffs; :meth:`resolve_case` still raises for a caller
            that asked about the case alone.
        """
        if submission_uploaded_date is None:
            # set arbitrary date if not provided
            submission_uploaded_date = metadata.submission.submission_date

            # Add submission date to ignore fields
            ignore_fields = ignore_fields or set()
            ignore_fields.add("submission_uploaded_date")

        # One transaction for the whole change set: the three parts describe one submission at
        # one moment, and computing them against separate snapshots let them disagree. Reading
        # the row here also means resolution and the field diff share it. ``get_submission``
        # would eagerly load the state log, which nothing in a diff reads.
        with self.transaction() as session:
            current_submission = session.get(Submission, submission_id)
            if current_submission is None:
                raise SubmissionNotFoundError(submission_id)

            try:
                case_link, case_link_error = (
                    self._diff_case_link(session, current_submission, metadata, ignore_fields or set()),
                    None,
                )
            except AmbiguousCaseError as exc:
                case_link, case_link_error = None, exc

            return SubmissionChangeSet(
                fields=self._diff_metadata(current_submission, metadata, submission_uploaded_date, ignore_fields),
                donors=self._diff_donors(session, submission_id, metadata),
                case_link=case_link,
                case_link_error=case_link_error,
            )

    def commit_changes(self, submission_id: str, changes: SubmissionChangeSet) -> None:
        """Write all pending changes of a change set to the database, as one transaction.

        A change set is one answer to one metadata file, so it is applied whole or not at all.
        Writing it piecemeal left a rejected write with the earlier parts already stored, and
        the operator was told only about the failure: re-running then previewed a different
        starting state than the one just shown.

        A pending case link in :attr:`SubmissionChangeSet.case_link` is applied via
        :meth:`assign_case`, creating the case first if none exists yet. It goes first, so
        that a case another process created in the meantime can be joined without discarding
        work already done here.

        :param submission_id: ID of the submission being updated.
        :param changes: Change set from :meth:`diff`.
        """
        with self.transaction() as session:
            if (link := changes.case_link) is not None:
                self.assign_case(
                    submission_id,
                    submitter_id=link.submitter_id,
                    local_case_id=link.local_case_id,
                    submission_type=link.submission_type,
                    session=session,
                )
            for field_diff in changes.fields.pending:
                self.modify_submission(submission_id, field_diff.key, field_diff.diff.after, session=session)
            # added/updated diffs carry a non-None `after`, deleted a non-None `before`,
            # by construction in Diff.classify
            for donor_diff in changes.donors.added:
                self.add_donor(cast(Donor, donor_diff.after), session=session)
            for donor_diff in changes.donors.updated:
                self.update_donor(cast(Donor, donor_diff.after), session=session)
            for donor_diff in changes.donors.deleted:
                self.delete_donor(cast(Donor, donor_diff.before), session=session)

    @staticmethod
    def _log_pending_changes(submission_id: str, changes: SubmissionChangeSet) -> None:
        """Emit info-level log lines summarising what is about to be committed."""
        sid = f"Submission: {submission_id}"

        pending_keys = [d.key for d in changes.fields.pending]
        unchanged_keys = [d.key for d in changes.fields.unchanged]
        if pending_keys:
            logger.info("%s - Updating fields: %s in database", sid, ", ".join(f'"{k}"' for k in pending_keys))
        if unchanged_keys:
            logger.info("%s - Not updating fields: %s in database", sid, ", ".join(f'"{k}"' for k in unchanged_keys))

        if (link := changes.case_link) is not None:
            target = f"existing case {link.after}" if link.after is not None else "a new case"
            logger.info(
                "%s - Linking to %s (submitter '%s', local case '%s')",
                sid,
                target,
                link.submitter_id,
                link.local_case_id,
            )

        donors = changes.donors
        if donors.unchanged:
            logger.info("%s - Keep existing donor(s): %s", sid, ", ".join(f'"{d.pseudonym}"' for d in donors.unchanged))
        if donors.added:
            logger.info(
                "%s - Adding new donor(s): %s",
                sid,
                ", ".join(f'"{d.pseudonym}"' for d in donors.added),
            )
        if donors.updated:
            logger.info(
                "%s - Modifying existing donor(s): %s",
                sid,
                ", ".join(f'"{d.pseudonym}"' for d in donors.updated),
            )
        if donors.deleted:
            logger.info("%s - Dropping donor(s): %s", sid, ", ".join(f'"{d.pseudonym}"' for d in donors.deleted))

    @staticmethod
    def assert_metadata_not_redacted(
        metadata: GrzSubmissionMetadata,
        submission_id: str,
        ignore_fields: set[str] | None = None,
    ) -> None:
        """Raise ``ValueError`` if ``metadata`` has redacted/missing required fields.

        Checks ``tan_g`` against :data:`REDACTED_TAN` and ``local_case_id``
        against :func:`is_redacted_local_case_id`. Each check can be
        bypassed by including the corresponding key in ``ignore_fields``:
        ``"tan_g"`` bypasses the redacted-TAN check; ``"pseudonym"`` bypasses
        the missing/redacted-``local_case_id`` check.

        :param metadata: Parsed submission metadata.
        :param submission_id: ID of the submission being populated.
        :param ignore_fields: Field keys whose check should be skipped.
        :raises ValueError: If ``tan_g`` or ``local_case_id`` is redacted/missing
            and the corresponding key is not in ``ignore_fields``.
        """
        ignore_fields = ignore_fields or set()
        if metadata.submission.tan_g == REDACTED_TAN and "tan_g" not in ignore_fields:
            raise ValueError(f"Submission {submission_id} has redacted tan_g in metadata.")
        if is_redacted_local_case_id(metadata.submission.local_case_id) and "pseudonym" not in ignore_fields:
            raise ValueError(f"Submission {submission_id} has missing or redacted local_case_id in metadata.")

    def populate(  # noqa: PLR0913
        self,
        submission_id: str,
        metadata: GrzSubmissionMetadata,
        submission_date: datetime.date | None,
        *,
        force: bool = False,
        on_missing: Literal["create", "error"] = "error",
        ignore_fields: set[str] | None = None,
    ) -> None:
        """Reconcile DB state for ``submission_id`` with ``metadata``.

        Rejects redacted ``tan_g`` or missing/redacted ``local_case_id`` via
        :meth:`assert_metadata_not_redacted` unless the corresponding key
        (``"tan_g"`` or ``"pseudonym"``) is in ``ignore_fields``. Computes diffs
        via :meth:`diff`, rejects destructive changes unless ``force``, and
        commits via :meth:`commit_changes`. Operational progress is logged via
        the module-level logger; callers configure verbosity through
        ``logging.getLogger("grz_db.models.submission")``.

        :param submission_id: ID of the submission being populated.
        :param metadata: Parsed submission metadata.
        :param submission_date: S3 last-modified date of the metadata file, if known.
        :param force: If ``True``, allow destructive updates/deletes of existing fields.
        :param on_missing: What to do when ``submission_id`` is not yet in the
            database. ``"error"`` (default) raises :class:`SubmissionNotFoundError`.
            ``"create"`` calls :meth:`add_submission` first and then proceeds; the
            resulting diff against the fresh row contains only additive changes,
            so ``force`` does not need to be set.
        :param ignore_fields: Field keys to skip during diff and redaction
            validation. See :meth:`assert_metadata_not_redacted`.
        :raises SubmissionNotFoundError: if the submission row is absent and
            ``on_missing`` is ``"error"``.
        :raises ValueError: if ``tan_g`` or ``local_case_id`` is redacted/missing
            and the corresponding key is not in ``ignore_fields``.
        :raises RuntimeError: if pending changes are destructive and ``force`` is False.
        """
        ignore_fields = ignore_fields or set()

        if not self.get_submission(submission_id):
            if on_missing == "error":
                raise SubmissionNotFoundError(submission_id)
            logger.warning("Submission %s does not exist. Creating ...", submission_id)
            self.add_submission(submission_id)
            logger.debug("Submission %s added to database.", submission_id)

        self.assert_metadata_not_redacted(metadata, submission_id, ignore_fields)

        changes = self.diff(submission_id, metadata, submission_date, ignore_fields=ignore_fields)

        if not force and changes.has_pending_destructive:
            raise RuntimeError(
                f"Would update/delete existing submission data "
                f"({', '.join(changes.destructive_changes)}) in the database, "
                f"but `force` not set. submission_id={submission_id!r}"
            )

        if changes.has_pending:
            logger.info("Submission: %s - Updating...", submission_id)
            self._log_pending_changes(submission_id, changes)
            self.commit_changes(submission_id, changes)
        else:
            logger.info("Submission: %s - No updates necessary.", submission_id)
