"""Shared fixtures for grz-db tests."""

import importlib.resources
import json

import pytest
from grz_db.testing import (  # noqa: F401
    _migrated_postgresql_template,
    _migrated_sqlite_template,
    _pg_test_dbnames,
    db,
    db_backend,
    db_test_connection,
    migrated_db_connection,
    postgresql_proc,
    test_author,
)
from grz_pydantic_models.submission.metadata import GrzSubmissionMetadata
from grz_pydantic_models_testing.example_metadata import grzctl as grzctl_metadata

TEST_METADATA_PATH = importlib.resources.files(grzctl_metadata).joinpath("metadata.json")


@pytest.fixture(scope="session")
def metadata() -> GrzSubmissionMetadata:
    """Load the wes_tumor_germline v1.2.1 example from grz-pydantic-models' own test resources."""
    with TEST_METADATA_PATH.open() as fh:
        return GrzSubmissionMetadata(**json.load(fh))
