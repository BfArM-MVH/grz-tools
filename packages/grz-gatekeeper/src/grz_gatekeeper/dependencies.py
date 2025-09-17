from fastapi import Depends
from grz_common.models.s3 import S3Options
from grz_db.models.author import Author
from grz_db.models.submission import SubmissionDb
from sqlmodel import Session

from .config import GatekeeperConfig, get_gatekeeper_config
from .session_store import gatekeeper_db_engine


def get_gatekeeper_db_session():
    """FastAPI dependency to get a session for the upload session store."""
    with Session(gatekeeper_db_engine) as session:
        yield session


def get_submission_db(config: GatekeeperConfig = Depends(get_gatekeeper_config)) -> SubmissionDb:
    """FastAPI dependency to get a SubmissionDb instance for business logic."""
    if path := config.db.author.private_key_path:
        with open(path, "rb") as f:
            private_key_bytes = f.read()
    elif key := config.db.author.private_key:
        private_key_bytes = key.encode("utf-8")
    else:
        raise ValueError("Either private_key or private_key_path must be provided.")

    author = Author(
        name=config.db.author.name,
        private_key_bytes=private_key_bytes,
        private_key_passphrase=config.db.author.private_key_passphrase,
    )
    return SubmissionDb(db_url=config.db.database_url, author=author)


def get_s3_options(config: GatekeeperConfig = Depends(get_gatekeeper_config)) -> S3Options:
    """FastAPI dependency to get S3 options from the configuration."""
    return config.s3
