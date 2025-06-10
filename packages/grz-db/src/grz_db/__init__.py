import datetime
import enum
import logging
import os
from collections.abc import Generator
from contextlib import contextmanager
from typing import Any, ClassVar, Generic, TypeVar

import cryptography
from alembic import command as alembic_command
from alembic.config import Config as AlembicConfig
from cryptography.hazmat.primitives.asymmetric.types import PrivateKeyTypes, PublicKeyTypes
from pydantic import ConfigDict
from sqlalchemy import JSON, Column
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm import selectinload
from sqlmodel import DateTime, Field, Relationship, Session, SQLModel, create_engine, select

__version__ = "0.1.0"

log = logging.getLogger(__name__)


class CaseInsensitiveStrEnum(enum.StrEnum):
    """
    A StrEnum that is case-insensitive for member lookup and comparison with strings.
    """

    @classmethod
    def _missing_(cls, value):
        """
        Override to allow case-insensitive lookup of enum members by value.
        e.g., MyEnum('value') will match MyEnum.VALUE.
        """
        if isinstance(value, str):
            for member in cls:
                if member.value.casefold() == value.casefold():
                    return member
        return None

    def __eq__(self, other):
        """
        Override to allow case-insensitive comparison of enum members by value.
        """
        if isinstance(other, enum.Enum):
            return self is other
        if isinstance(other, str):
            return self.value.casefold() == other.casefold()
        return NotImplemented

    def __hash__(self):
        """
        Override to make hash consistent with eq.
        """
        return hash(self.value.casefold())


def serialize_datetime_to_iso_z(dt: datetime.datetime) -> str:
    """
    Serializes a datetime object to a canonical ISO 8601 string format with 'Z' for UTC.
    """
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=datetime.UTC)

    if dt.tzinfo != datetime.UTC and dt.utcoffset() != datetime.timedelta(0):
        dt = dt.astimezone(datetime.UTC)

    return dt.isoformat()


class ListableEnum(enum.StrEnum):
    """Mixin for enum classes whose members can be listed."""

    @classmethod
    def list(cls) -> list[str]:
        """Returns a list of enum members."""
        return list(map(lambda c: c.value, cls))


class SubmissionStateEnum(CaseInsensitiveStrEnum, ListableEnum):
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


class ChangeRequestEnum(CaseInsensitiveStrEnum, ListableEnum):
    """Change request enum."""

    MODIFY = "Modify"
    DELETE = "Delete"
    TRANSFER = "Transfer"


class BaseSignablePayload(SQLModel):
    """
    Base class for SQLModel based payloads
    that can be signed and can be converted to bytes for verification.
    Provides a default `to_bytes` method using pydantic's JSON serialization.
    Provides a default `sign` method using the private key of the author.
    """

    model_config = ConfigDict(
        json_encoders={datetime.datetime: serialize_datetime_to_iso_z},
        populate_by_name=True,
    )

    def to_bytes(self) -> bytes:
        """
        Default serialization: JSON string encoded to UTF-8.
        """
        payload_json = self.model_dump_json(by_alias=True)
        return payload_json.encode("utf8")

    def sign(self, private_key: PrivateKeyTypes) -> bytes:
        """Sign this payload using the given private key."""
        bytes_to_sign = self.to_bytes()
        signature = private_key.sign(bytes_to_sign)
        public_key_of_private = private_key.public_key()
        public_key_of_private.verify(signature, bytes_to_sign)
        return signature


P = TypeVar("P", bound=BaseSignablePayload)


class VerifiableLog(Generic[P]):
    """
    Mixin class for SQLModels that store a signature and can be verified.
    Subclasses MUST:
    1. Define `payload_model_class: ClassVar[type[P]]`.
    2. Have an instance attribute `signature: str`.
    """

    signature: str
    payload_model_class: ClassVar[type[P]]

    def __init_subclass__(cls, **kwargs: Any) -> None:
        super().__init_subclass__(**kwargs)
        if not hasattr(cls, "payload_model_class"):
            raise TypeError(f"Class {cls.__name__} lacks 'payload_model_class' attribute required by VerifiableLog.")
        if not (isinstance(cls.payload_model_class, type) and issubclass(cls.payload_model_class, BaseSignablePayload)):
            raise TypeError(
                f"'payload_model_class' in {cls.__name__} must be a class and a subclass of BaseSignedPayload. "
                f"Got: {cls.payload_model_class}"
            )

    def verify(self, public_key: PublicKeyTypes) -> bool:
        """Verify the signature of this log entry."""
        if not hasattr(self, "signature") or not isinstance(self.signature, str) or not self.signature:
            log.warning(f"Missing/invalid signature for {self.__class__.__name__} (id: {getattr(self, 'id', 'N/A')}).")
            return False

        signature_bytes = bytes.fromhex(self.signature)
        data_for_payload = self.model_dump(by_alias=True, exclude={"signature", "payload_model_class"})
        payload_to_verify = self.payload_model_class(**data_for_payload)
        bytes_to_verify = payload_to_verify.to_bytes()

        try:
            public_key.verify(signature_bytes, bytes_to_verify)
        except cryptography.exceptions.InvalidSignature:
            return False
        except:
            raise
        return True


class SubmissionBase(SQLModel):
    """Submission base model."""

    tan_g: str | None = Field(default=None, unique=True, index=True, alias="tanG")
    pseudonym: str | None = Field(default=None, index=True)


class Submission(SubmissionBase, table=True):
    """Submission table model."""

    __tablename__ = "submissions"

    id: str = Field(primary_key=True, index=True)

    states: list["SubmissionStateLog"] = Relationship(back_populates="submission")

    changes: list["ChangeRequestLog"] = Relationship(back_populates="submission")


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

    model_config = ConfigDict(
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

    payload_model_class = SubmissionStateLogPayload

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

    model_config = ConfigDict(
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

    payload_model_class = ChangeRequestLogPayload

    id: int | None = Field(default=None, primary_key=True, index=True)
    submission_id: str = Field(foreign_key="submissions.id", index=True)

    author_name: str = Field(index=True)
    signature: str

    submission: Submission | None = Relationship(back_populates="changes")


class SubmissionNotFoundError(ValueError):
    """Exception for when a submission is not found in the database."""

    def __init__(self, submission_id: str):
        super().__init__(f"Submission not found for ID {submission_id}")


class DuplicateSubmissionError(ValueError):
    """Exception for when a submission ID already exists in the database."""

    def __init__(self, submission_id: str):
        super().__init__(f"Duplicate submission ID {submission_id}")


class DuplicateTanGError(ValueError):
    """Exception for when a tanG is already in use."""

    def __init__(self, tan_g: str):
        super().__init__(f"Duplicate tanG {tan_g}")


class DatabaseConfigurationError(Exception):
    """Exception for database configuration issues."""

    pass


class Author:
    def __init__(self, name: str, private_key_bytes: bytes):
        self.name = name
        self.private_key_bytes = private_key_bytes

    def private_key(self) -> PrivateKeyTypes:
        from functools import partial
        from getpass import getpass

        from cryptography.hazmat.primitives.serialization import load_ssh_private_key

        passphrase = os.getenv("GRZ_DB_AUTHOR_PASSPHRASE")
        passphrase_callback = (lambda: passphrase) if passphrase else None

        if not passphrase:
            passphrase_callback = partial(getpass, prompt=f"Passphrase for GRZ DB author ({self.name}'s) private key: ")
        log.info(f"Loading private key of {self.name}…")
        private_key = load_ssh_private_key(
            self.private_key_bytes,
            password=passphrase_callback().encode("utf-8"),
        )
        return private_key


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
    def get_session(self) -> Generator[Session, Any, None]:
        """Get an sqlmodel session."""
        with Session(self.engine) as session:
            yield session

    def _get_alembic_config(self, alembic_ini_path: str) -> AlembicConfig:
        """
        Loads the alembic configuration.

        Args:
            alembic_ini_path: Path to alembic ini file.
        """
        if not alembic_ini_path or not os.path.exists(alembic_ini_path):
            raise ValueError(f"Alembic configuration file not found at: {alembic_ini_path}")

        alembic_cfg = AlembicConfig(alembic_ini_path)
        alembic_cfg.set_main_option("sqlalchemy.url", str(self.engine.url))
        alembic_cfg.set_main_option("script_location", "grz_db:migrations")
        return alembic_cfg

    def initialize_schema(self):
        """Initialize the database."""
        SQLModel.metadata.create_all(self.engine, checkfirst=True)

    def upgrade_schema(self, alembic_ini_path: str, revision: str = "head"):
        """
        Upgrades the database schema using alembic.

        Args:
            alembic_ini_path: Path to the alembic.ini file.
            revision: The Alembic revision to upgrade to (default: 'head').

        Raises:
            RuntimeError: For underlying Alembic errors.
        """
        alembic_cfg = self._get_alembic_config(alembic_ini_path)
        try:
            alembic_command.upgrade(alembic_cfg, revision)
        except Exception as e:
            raise RuntimeError(f"Alembic upgrade failed: {e}") from e

    def add_submission(
        self,
        submission_id: str,
        tan_g: str | None = None,
        pseudonym: str | None = None,
    ) -> Submission:
        """
        Adds a submission to the database.

        Args:
            submission_id: Submission ID.
            tan_g: tanG if in phase 0
            pseudonym: pseudonym if phase >= 0

        Returns:
            An instance of Submission.
        """
        with self.get_session() as session:
            existing_submission = session.get(Submission, submission_id)
            if existing_submission:
                raise DuplicateSubmissionError(submission_id)

            submission_create = SubmissionCreate(id=submission_id, tan_g=tan_g, pseudonym=pseudonym)
            db_submission = Submission.model_validate(submission_create)

            session.add(db_submission)
            try:
                session.commit()
                session.refresh(db_submission)
                return db_submission
            except IntegrityError as e:
                session.rollback()
                if "UNIQUE constraint failed: submissions.tanG" in str(e) and tan_g:
                    raise DuplicateTanGError(tan_g) from e
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
        with self.get_session() as session:
            submission = session.get(Submission, submission_id)
            if not submission:
                raise SubmissionNotFoundError(submission_id)

            state_log_payload = SubmissionStateLogPayload(
                submission_id=submission_id, author_name=self._author.name, state=state, data=data
            )
            signature = state_log_payload.sign(self._author.private_key())

            state_log_create = SubmissionStateLogCreate(**state_log_payload.model_dump(), signature=signature.hex())
            db_state_log = SubmissionStateLog.model_validate(state_log_create)
            session.add(db_state_log)

            # Remove tanG once it has been reported?
            if state == SubmissionStateEnum.REPORTED and submission.tan_g is not None:
                submission.tan_g = None
                session.add(submission)

            try:
                session.commit()
                session.refresh(db_state_log)
                if state == SubmissionStateEnum.REPORTED and submission.tan_g is None:
                    session.refresh(submission)
                return db_state_log
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
        with self.get_session() as session:
            statement = (
                select(Submission).where(Submission.id == submission_id).options(selectinload(Submission.states))
            )
            submission = session.exec(statement).first()
            return submission

    def list_submissions(self) -> list[Submission]:
        """
        Lists all submissions in the database.

        Returns:
            A list of all submissions in the database, ordered by their ID.
        """
        with self.get_session() as session:
            statement = select(Submission).options(selectinload(Submission.states)).order_by(Submission.id)
            submissions = session.exec(statement).all()
            return submissions
