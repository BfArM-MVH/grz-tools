import datetime
from collections.abc import Generator, Sequence
from contextlib import contextmanager
from operator import attrgetter
from typing import Any, ClassVar, Optional

import sqlalchemy as sa
from alembic import command as alembic_command
from alembic.config import Config as AlembicConfig
from alembic.runtime.migration import MigrationContext
from alembic.script import ScriptDirectory as AlembicScriptDirectory
from grz_pydantic_models.submission.metadata import (
    CoverageType,
    DiseaseType,
    GenomicDataCenterId,
    GenomicStudySubtype,
    GenomicStudyType,
    LibraryType,
    Relation,
    ResearchConsentNoScopeJustification,
    SequenceSubtype,
    SequenceType,
    SubmissionType,
    SubmitterId,
    Tan,
)
from pydantic import ConfigDict, field_serializer, field_validator
from sqlalchemy import JSON, Column
from sqlalchemy import func as sqlfn
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm import selectinload
from sqlmodel import DateTime, Field, Relationship, Session, SQLModel, create_engine, select

from ..common import (
    CaseInsensitiveStrEnum,
    ListableEnum,
    serialize_datetime_to_iso_z,
)
from ..errors import DuplicateSubmissionError, DuplicateTanGError, SubmissionNotFoundError
from .author import Author
from .base import BaseSignablePayload, VerifiableLog


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
    REPORTED = "Reported"
    QCING = "QCing"
    QCED = "QCed"
    CLEANING = "Cleaning"
    CLEANED = "Cleaned"
    FINISHED = "Finished"
    ERROR = "Error"


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
    submission_date: datetime.date | None = None
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


class Submission(SubmissionBase, table=True):
    """Submission table model."""

    __tablename__ = "submissions"
    __table_args__ = {"extend_existing": True}

    id: str = Field(primary_key=True, index=True)

    @field_validator("id")
    @classmethod
    def validate_id_pattern(cls, v: str) -> str:
        import re

        pattern = r"^[0-9]{9}_\d{4}-\d{2}-\d{2}_[a-f0-9]{8}$"
        if not re.match(pattern, v):
            raise ValueError(f"Submission ID '{v}' does not match the required pattern.")
        return v

    states: list["SubmissionStateLog"] = Relationship(back_populates="submission")

    changes: list["ChangeRequestLog"] = Relationship(back_populates="submission")

    def get_latest_state(self, filter_to_type: SubmissionStateEnum | None = None) -> Optional["SubmissionStateLog"]:
        states = filter(lambda state: state.state == filter_to_type, self.states) if filter_to_type else self.states
        states = sorted(states, key=attrgetter("timestamp"))
        return states[-1] if states else None


class SubmissionStateLogBase(SQLModel):
    """
    Submission state log base model.
    Holds state information for each submission.
    Timestamped.
    Can optionally have associated JSON data.
    """

    state: SubmissionStateEnum
    data: dict[str, Any] | None = Field(default=None, sa_column=Column(JSON))
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


class ChangeRequestLogBase(SQLModel):
    """
    Base model for change request logs.
    Timestamped.
    Can optionally have associated JSON data.
    """

    change: ChangeRequestEnum
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


class ChangeRequestLogPayload(ChangeRequestLogBase, BaseSignablePayload):
    """
    Used to bundle data for signature calculation.
    """

    submission_id: str
    author_name: str


class ChangeRequestLog(ChangeRequestLogBase, VerifiableLog[ChangeRequestLogPayload], table=True):
    """Change-request log table model."""

    __tablename__ = "submission_change_requests"
    __table_args__ = {"extend_existing": True}

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
    relation: Relation
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

    model_config = ConfigDict(  # type: ignore
        populate_by_name=True,
    )

    @field_serializer("timestamp")
    def serialize_timestamp(self, ts: datetime.datetime) -> str:
        return serialize_datetime_to_iso_z(ts)


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

    @contextmanager
    def _get_session(self) -> Generator[Session, Any, None]:
        """Get an sqlmodel session."""
        if not self._at_latest_schema():
            raise OutdatedDatabaseSchemaError(
                "Database not at latest schema. Please backup the database and then attempt a migration with `grzctl db upgrade`."
            )
        with Session(self.engine) as session:
            yield session

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
        with self._get_session() as session:
            existing_submission = session.get(Submission, submission_id)
            if existing_submission:
                raise DuplicateSubmissionError(submission_id)

            submission_create = SubmissionCreate(id=submission_id)
            db_submission = Submission.model_validate(submission_create)

            session.add(db_submission)
            try:
                session.commit()
                session.refresh(db_submission)
                return db_submission
            except IntegrityError as e:
                session.rollback()
                raise e
            except Exception:
                session.rollback()
                raise

    def modify_submission(self, submission_id: str, key: str, value: str) -> Submission:
        if key not in SubmissionBase.model_fields:
            raise ValueError(f"Unknown column key '{key}'")
        elif key in SubmissionBase.immutable_fields:
            raise ValueError(f"Column '{key}' is read-only and cannot be modified.")

        with self._get_session() as session:
            submission = session.get(Submission, submission_id)
            if submission is None:
                raise SubmissionNotFoundError(submission_id)

            setattr(submission, key, value)
            session.add(submission)
            try:
                session.commit()
                session.refresh(submission)
                return submission
            except IntegrityError as e:
                session.rollback()
                if "UNIQUE constraint failed: submissions.tanG" in str(e) and key == "tan_g":
                    raise DuplicateTanGError() from e
                raise
            except Exception:
                session.rollback()
                raise

    def update_submission_state(
        self,
        submission_id: str,
        state: SubmissionStateEnum,
        data: dict | None = None,
    ) -> SubmissionStateLog:
        """
        Updates a submission's state to the specified state.

        Args:
            submission_id: Submission ID of the submission to update.
            state: New state of the submission.
            data: Optional data to attach to the update.

        Returns:
            An instance of SubmissionStateLog.
        """
        with self._get_session() as session:
            submission = session.get(Submission, submission_id)
            if not submission:
                raise SubmissionNotFoundError(submission_id)
            if not self._author:
                raise ValueError("No author defined")

            state_log_payload = SubmissionStateLogPayload(
                submission_id=submission_id, author_name=self._author.name, state=state, data=data
            )
            signature = state_log_payload.sign(self._author.private_key())

            state_log_create = SubmissionStateLogCreate(**state_log_payload.model_dump(), signature=signature.hex())
            db_state_log = SubmissionStateLog.model_validate(state_log_create)
            session.add(db_state_log)

            try:
                session.commit()
                session.refresh(db_state_log)
                return db_state_log
            except Exception:
                session.rollback()
                raise

    def get_donors(self, submission_id: str, pseudonym: str | None = None) -> tuple[Donor, ...]:
        """Retrieve all donors for a given submission, or, optionally, only for a specific pseudonym."""
        with self._get_session() as session:
            statement = select(Donor).where(Donor.submission_id == submission_id)
            if pseudonym is not None:
                statement = statement.where(Donor.pseudonym == pseudonym)
            donors = tuple(session.exec(statement).all())
        return donors

    def add_donor(self, donor: Donor) -> Donor:
        """Add a donor to the database."""
        with self._get_session() as session:
            session.add(donor)

            try:
                session.commit()
                session.refresh(donor)
                return donor
            except Exception as e:
                session.rollback()
                raise e

    def update_donor(self, updated_donor: Donor) -> Donor:
        """Update a donor in the database."""
        with self._get_session() as session:
            statement = (
                select(Donor)
                .where(Donor.submission_id == updated_donor.submission_id)
                .where(Donor.pseudonym == updated_donor.pseudonym)
            )
            db_donor = session.exec(statement).first()

            if db_donor is None:
                raise RuntimeError("Cannot update a donor that doesn't yet exist in the database.")

            if db_donor == updated_donor:
                # nothing to do
                return db_donor

            for field in db_donor.model_fields:
                old_value = getattr(db_donor, field)
                new_value = getattr(updated_donor, field)
                if old_value != new_value:
                    setattr(db_donor, field, new_value)

            session.add(db_donor)

            try:
                session.commit()
                session.refresh(db_donor)
                return db_donor
            except Exception as e:
                session.rollback()
                raise e

    def delete_donor(self, donor: Donor) -> None:
        """Delete a donor from the database."""
        with self._get_session() as session:
            session.delete(donor)

            try:
                session.commit()
            except Exception as e:
                session.rollback()
                raise e

    def get_detailed_qc_results(self, submission_id: str) -> tuple[DetailedQCResult, ...]:
        """Retrieve all detailed QC results for a given submission."""
        with self._get_session() as session:
            statement = select(DetailedQCResult).where(DetailedQCResult.submission_id == submission_id)
            results = tuple(session.exec(statement).all())
        return results

    def add_detailed_qc_result(self, result: DetailedQCResult) -> DetailedQCResult:
        """Add or update a detailed QC result to/in the database."""
        with self._get_session() as session:
            session.add(result)

            try:
                session.commit()
                session.refresh(result)
                return result
            except Exception as e:
                session.rollback()
                raise e

    def add_change_request(
        self,
        submission_id: str,
        change: ChangeRequestEnum,
        data: dict | None = None,
    ) -> ChangeRequestLog:
        """
        Register a change request for a submission.

        Args:
            submission_id: Submission ID of the submission to register a change request for.
            change: Requested change.
            data: Optional data to attach to the update.

        Returns:
            An instance of ChangeRequestLog.
        """
        with self._get_session() as session:
            submission = session.get(Submission, submission_id)
            if not submission:
                raise SubmissionNotFoundError(submission_id)
            if not self._author:
                raise ValueError("No author defined")

            change_request_log_payload = ChangeRequestLogPayload(
                submission_id=submission_id, author_name=self._author.name, change=change, data=data
            )
            signature = change_request_log_payload.sign(self._author.private_key())

            change_request_log_create = ChangeRequestLogCreate(
                **change_request_log_payload.model_dump(), signature=signature.hex()
            )
            db_change_request_log = ChangeRequestLog.model_validate(change_request_log_create)
            session.add(db_change_request_log)

            try:
                session.commit()
                session.refresh(db_change_request_log)
                return db_change_request_log
            except Exception:
                session.rollback()
                raise

    def get_submission(self, submission_id: str) -> Submission | None:
        """
        Retrieves a submission and its state history.

        Args:
            submission_id: Submission ID of the submission to retrieve.

        Returns:
            An instance of Submission or None.
        """
        with self._get_session() as session:
            statement = (
                select(Submission).where(Submission.id == submission_id).options(selectinload(Submission.states))  # type: ignore[arg-type]
            )
            submission = session.exec(statement).first()
            return submission

    def list_submissions(self, limit: int | None) -> Sequence[Submission]:
        """
        Lists all submissions in the database.

        Returns:
            A list of all submissions in the database. Ordered by latest
            submission state timestamp if not null, otherwise use submission
            date, with submissions missing both of these sorting first.
        """
        with self._get_session() as session:
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
                    sqlfn.coalesce(latest_state_per_submission.c.timestamp, Submission.submission_date)
                    .desc()
                    .nulls_first()
                )
            )
            if limit is not None:
                statement = statement.limit(limit)
            submissions = session.exec(statement).all()
            return submissions

    def list_processed_between(self, start: datetime.date, end: datetime.date) -> Sequence[Submission]:
        """
        Lists all submissions processed between the given start and end dates, inclusive.
        Processed is defined as either reported (Prüfbericht submitted) or detailed QC finished.
        """
        with self._get_session() as session:
            reported_within_window = (
                select(SubmissionStateLog.submission_id)
                .where(SubmissionStateLog.state.in_([SubmissionStateEnum.REPORTED, SubmissionStateEnum.QCED]))  # type: ignore[attr-defined]
                .where(SubmissionStateLog.timestamp.between(start, end))  # type: ignore[attr-defined]
                .subquery()
            )
            statement = (
                select(Submission)
                .options(selectinload(Submission.states))  # type: ignore[arg-type]
                .join(reported_within_window, Submission.id == reported_within_window.c.submission_id)  # type: ignore[arg-type]
                .distinct()
            )
            submissions = session.exec(statement).all()
            return submissions

    def list_change_requests(self) -> Sequence[Submission]:
        """
        Lists all submissions in the database.

        Returns:
            A list of all submissions in the database, ordered by their ID.
        """
        with self._get_session() as session:
            statement = (
                select(Submission)
                .where(Submission.changes.any())  # type: ignore[attr-defined]
                .options(selectinload(Submission.changes))  # type: ignore[arg-type]
                .order_by(Submission.id)
            )
            change_requests = session.exec(statement).all()
            return change_requests
