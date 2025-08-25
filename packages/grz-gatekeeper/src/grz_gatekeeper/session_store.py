from datetime import UTC, datetime

from grz_pydantic_models.gatekeeper import FileUploadInfo
from sqlalchemy import JSON, Column
from sqlmodel import Field, SQLModel, create_engine

GATEKEEPER_DATABASE_URL = "sqlite:///gatekeeper.sessions.db"  # TODO: expose in config, use ":memory:" as default


class UploadSession(FileUploadInfo, SQLModel, table=True):
    id: int | None = Field(default=None, primary_key=True)
    submission_id: str = Field(index=True)
    created_at: datetime = Field(default_factory=lambda: datetime.now(UTC))
    last_modified_at: datetime = Field(default_factory=lambda: datetime.now(UTC))
    parts: list = Field(sa_column=Column(JSON))


gatekeeper_db_engine = create_engine(GATEKEEPER_DATABASE_URL)


def init_gatekeeper_sessions_db():
    SQLModel.metadata.create_all(gatekeeper_db_engine)
