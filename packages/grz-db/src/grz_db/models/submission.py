import datetime
from collections.abc import Generator, Sequence
from contextlib import contextmanager
from operator import attrgetter
from typing import Annotated, Any, ClassVar, Optional

from alembic import command as alembic_command
from alembic.config import Config as AlembicConfig
from alembic.runtime.migration import MigrationContext
from alembic.script import ScriptDirectory as AlembicScriptDirectory
from grz_pydantic_models.submission.metadata import (
    DiseaseType,
    GenomicDataCenterId,
    LibraryType,
    SubmissionType,
    SubmitterId,
    Tan,
)
from pydantic import ConfigDict, StringConstraints
from sqlalchemy import JSON, Column
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
    disease_type: DiseaseType | None = None
    library_type: LibraryType | None = None
    basic_qc_passed: bool | None = None

    # fields also for Tätigkeitsbericht
    consented: bool | None = None
    detailed_qc_passed: bool | None = None


class Submission(SubmissionBase, table=True):
    """Submission table model."""

    __tablename__ = "submissions"

    id: Annotated[str, StringConstraints(pattern=r"^[0-9]{9}_\d{4}-\d{2}-\d{2}_[a-f0-9]{8}$")] = Field(
        primary_key=True, index=True
    )

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
        json_encoders={datetime.datetime: serialize_datetime_to_iso_z},
        populate_by_name=True,
    )


class SubmissionStateLogPayload(SubmissionStateLogBase, BaseSignablePayload):
    """
    Used to bundle data for signature calculation.
    """

    submission_id: str
    author_name: str


class SubmissionStateLog(SubmissionStateLogBase, VerifiableLog[SubmissionStateLogPayload], table=True):
    """Submission state log table model."""

    __tablename__ = "submission_states"

    _payload_model_class = SubmissionStateLogPayload

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
        json_encoders={datetime.datetime: serialize_datetime_to_iso_z},
        populate_by_name=True,
    )


class ChangeRequestLogPayload(ChangeRequestLogBase, BaseSignablePayload):
    """
    Used to bundle data for signature calculation.
    """

    submission_id: str
    author_name: str


class ChangeRequestLog(ChangeRequestLogBase, VerifiableLog[ChangeRequestLogPayload], table=True):
    """Change-request log table model."""

    __tablename__ = "submission_change_requests"

    _payload_model_class = ChangeRequestLogPayload

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
            A list of all submissions in the database, ordered by their submission date, latest first.
        """
        with self._get_session() as session:
            statement = (
                select(Submission).options(selectinload(Submission.states)).order_by(Submission.submission_date.desc())  # type: ignore[arg-type,union-attr]
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
