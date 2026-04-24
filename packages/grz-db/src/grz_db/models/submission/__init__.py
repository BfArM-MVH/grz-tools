import calendar
import datetime
import logging
import math
import random
import re
from collections.abc import Generator, Sequence
from contextlib import contextmanager
from operator import attrgetter
from typing import Any, ClassVar, Optional, Self

import sqlalchemy as sa
import sqlalchemy.dialects.postgresql as sa_psql
from alembic import command as alembic_command
from alembic.config import Config as AlembicConfig
from alembic.runtime.migration import MigrationContext
from alembic.script import ScriptDirectory as AlembicScriptDirectory
from grz_pydantic_models.dates import date_to_quarter_year, quarter_date_bounds
from grz_pydantic_models.submission.metadata import (
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
)
from grz_pydantic_models.submission.metadata.v1 import Donor as MetadataDonor
from pydantic import ConfigDict, field_serializer, field_validator
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
    DuplicateSubmissionError,
    DuplicateTanGError,
    SubmissionBasicQCNotPassedError,
    SubmissionDateIsNoneError,
    SubmissionNotFoundError,
    SubmissionTypeIsNoneError,
)
from ..author import Author
from ..base import BaseSignablePayload, VerifiableLog
from .diff import (  # noqa: F401
    Diff,
    DiffState,
    DonorDiff,
    DonorsDiffCollection,
    FieldDiff,
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

    # Database column indicating whether a submission is selected for in-depth QC (True/False) or not yet decided (None).
    selected_for_qc: bool | None = None

    # extra fields
    submission_size: int | None = Field(default=None, sa_type=BigInteger)
    submission_metadata: dict[str, Any] | None = Field(
        default=None,
        sa_column=Column(JSON().with_variant(sa_psql.JSONB, "postgresql")),
    )


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

    states: list["SubmissionStateLog"] = Relationship(back_populates="submission")

    changes: list["ChangeRequestLog"] = Relationship(back_populates="submission")

    def get_latest_state(self, filter_to_type: SubmissionStateEnum | None = None) -> Optional["SubmissionStateLog"]:
        states = filter(lambda state: state.state == filter_to_type, self.states) if filter_to_type else self.states
        states = sorted(states, key=attrgetter("timestamp"))
        return states[-1] if states else None

    @classmethod
    def from_metadata(
        cls,
        submission_id: str,
        metadata: "GrzSubmissionMetadata",
        submission_date: datetime.date | None,
    ) -> "Submission":
        """Construct a Submission populated with values derived from parsed metadata.

        Only the fields that can be sourced from metadata are set; system-managed
        fields (e.g. ``basic_qc_passed``, ``selected_for_qc``) are left at their
        defaults so that ``model_fields_set`` reliably indicates which fields to
        compare during a diff.
        """
        return cls.model_validate(
            {
                "id": submission_id,
                "tan_g": metadata.submission.tan_g,
                "submission_type": metadata.submission.submission_type,
                "submitter_id": metadata.submission.submitter_id,
                "coverage_type": metadata.submission.coverage_type,
                "disease_type": metadata.submission.disease_type,
                "genomic_study_type": metadata.submission.genomic_study_type,
                "genomic_study_subtype": metadata.submission.genomic_study_subtype,
                "pseudonym": metadata.submission.local_case_id,
                "data_node_id": metadata.submission.genomic_data_center_id,
                "consented": metadata.consents_to_research(date=datetime.date.today()),
                "submission_size": metadata.get_submission_size(),
                "submission_date": submission_date
                if submission_date is not None
                else metadata.submission.submission_date,
                "submission_metadata": metadata.to_redacted_dict(),
            }
        )


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
    ) -> Self:
        return cls.model_validate(
            dict(
                submission_id=submission_id,
                pseudonym="index" if donor.relation == Relation.index_ else donor.donor_pseudonym,
                relation=Relation(donor.relation),
                library_types={datum.library_type for datum in donor.lab_data},
                sequence_types={datum.sequence_type for datum in donor.lab_data},
                sequence_subtypes={datum.sequence_subtype for datum in donor.lab_data},
                mv_consented=donor.consents_to_mv(),
                research_consented=donor.consents_to_research(date=datetime.date.today()),
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

    def modify_submission(self, submission_id: str, key: str, value: Any) -> Submission:  # noqa: C901
        if key not in SubmissionBase.model_fields:
            raise ValueError(f"Unknown column key '{key}'")
        elif key in SubmissionBase.immutable_fields:
            raise ValueError(f"Column '{key}' is read-only and cannot be modified.")

        with self._get_session() as session:
            submission = session.get(Submission, submission_id)
            if submission is None:
                raise SubmissionNotFoundError(submission_id)

            setattr(submission, key, value)
            if key == "basic_qc_passed":
                # Basic QC state changed -> Align the in-depth QC queue to the new state
                queue_entry = session.get(QCQueueEntry, submission_id)

                if submission.basic_qc_passed is True and queue_entry is None:
                    # Basic QC passed -> Ensure that submission is tracked in the in-depth QC queue
                    session.add(QCQueueEntry(submission_id=submission_id))
                elif submission.basic_qc_passed is not True and queue_entry is not None:
                    # Basic QC failed -> Ensure that submission is absent from the in-depth QC queue
                    session.delete(queue_entry)

                # Keep selection flag aligned with failed basic QC.
                if submission.basic_qc_passed is False:
                    submission.selected_for_qc = False
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

    def update_submission(self, submission: Submission) -> Submission:
        """
        Persists changes made to a Submission object back to the database.

        :param submission: The Submission instance with updated field values.
        :return: The updated Submission instance.
        """
        with self._get_session() as session:
            db_submission = session.get(Submission, submission.id)
            if db_submission is None:
                raise SubmissionNotFoundError(submission.id)

            for field in SubmissionBase.model_fields:
                if field in SubmissionBase.immutable_fields:
                    continue
                setattr(db_submission, field, getattr(submission, field))

            session.add(db_submission)
            try:
                session.commit()
                session.refresh(db_submission)
                return db_submission
            except IntegrityError as e:
                session.rollback()
                if "UNIQUE constraint failed: submissions.tanG" in str(e):
                    raise DuplicateTanGError() from e
                raise
            except Exception:
                session.rollback()
                raise

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
        with self._get_session() as session:
            return session.exec(
                select(Submission)
                .options(selectinload(Submission.states))  # type: ignore[arg-type]
                .join(QCQueueEntry, QCQueueEntry.submission_id == Submission.id)  # type: ignore[arg-type]
                .where(Submission.submission_type == SubmissionType.initial)
                .where(Submission.basic_qc_passed)  # type: ignore[arg-type]
                .where(Submission.submission_date.between(start_date, end_date))  # type: ignore[union-attr]
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
        submission_quarter, submission_year = date_to_quarter_year(submission.submission_date)  # type: ignore[arg-type]
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

            for field in Donor.model_fields:
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
        with self._get_session() as session:
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
        with self._get_session() as session:
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
        submission_date = submission.submission_date
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

    def _diff_metadata(
        self,
        submission_id: str,
        metadata: "GrzSubmissionMetadata",
        submission_date: datetime.date | None,
        ignore_fields: set[str] | None = None,
    ) -> SubmissionDiffCollection:
        """Compare a submission's current database state against fresh metadata.

        :param submission_id: Submission ID to look up.
        :param metadata: Parsed metadata from the submission's ``metadata.json``.
        :param submission_date: Explicit submission date; falls back to the value in *metadata* when ``None``.
        :param ignore_fields: Field names to skip entirely during the comparison.
        :returns: A :class:`SubmissionDiffCollection` instance summarising all detected differences.
        """
        db_submission = self.get_submission(submission_id)
        if db_submission is None:
            raise SubmissionNotFoundError(submission_id)

        if ignore_fields is None:
            ignore_fields = set()

        fresh_submission = Submission.from_metadata(submission_id, metadata, submission_date)

        result = SubmissionDiffCollection()
        for key in fresh_submission.model_fields_set - (ignore_fields or set()):
            field_diff = FieldDiff.classify_field(key, getattr(db_submission, key), getattr(fresh_submission, key))
            result.append(field_diff)

        return result

    def _diff_donors(
        self,
        submission_id: str,
        metadata: "GrzSubmissionMetadata",
    ) -> DonorsDiffCollection:
        """Diff all donors in *metadata* against the current database state.

        :param submission_id: Submission ID to look up donors for.
        :param metadata: Parsed metadata from the submission's ``metadata.json``.
        :returns: A fully populated :class:`DonorDiff`.
        """
        donors_in_db_submission = {donor.pseudonym: donor for donor in self.get_donors(submission_id=submission_id)}
        donors_in_metadata = {
            (d := Donor.from_donor_metadata(submission_id, donor)).pseudonym: d for donor in metadata.donors
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
        metadata: "GrzSubmissionMetadata",
        submission_date: datetime.date | None,
        ignore_fields: set[str] | None = None,
    ) -> tuple[SubmissionDiffCollection, DonorsDiffCollection]:
        submission_diff = self._diff_metadata(submission_id, metadata, submission_date, ignore_fields)
        donor_diff = self._diff_donors(submission_id, metadata)
        return submission_diff, donor_diff

    def commit_changes(
        self,
        submission_id: str,
        submission_diff: SubmissionDiffCollection | None = None,
        donors_diff: DonorsDiffCollection | None = None,
    ) -> None:
        """Write all pending metadata and donor diffs to the database.
        Can be obtained by calling :func:`diff`

        :param db: Database service instance to write to.
        :param submission_id: ID of the submission being updated.
        :param submission_diff: Diff result from :func:`diff_metadata`.
        :param donors_diff: Diff result from :func:`build_donor_diff`.
        """
        if submission_diff is not None:
            for field_diff in submission_diff.pending:
                self.modify_submission(submission_id, field_diff.key, field_diff.diff.after)
        if donors_diff is not None:
            for donor_diff in donors_diff.added:
                assert donor_diff.after is not None, "Added NoneType donor, this should not happen"  # noqa: S101
                self.add_donor(donor_diff.after)
            for donor_diff in donors_diff.updated:
                assert donor_diff.after is not None, "Updated NoneType donor, this should not happen"  # noqa: S101
                self.update_donor(donor_diff.after)
            for donor_diff in donors_diff.deleted:
                assert donor_diff.before is not None, "Removed NoneType donor, this should not happen"  # noqa: S101
                self.delete_donor(donor_diff.before)
