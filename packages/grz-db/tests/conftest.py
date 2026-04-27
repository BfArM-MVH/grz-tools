"""Shared fixtures for grz-db tests."""

import importlib.resources
import json
from shutil import which

import psycopg
import pytest
from grz_pydantic_models.submission.metadata import GrzSubmissionMetadata
from grz_pydantic_models_testing.example_metadata import grzctl as grzctl_metadata

TEST_METADATA_PATH = importlib.resources.files(grzctl_metadata).joinpath("metadata.json")


@pytest.fixture(scope="session")
def metadata() -> GrzSubmissionMetadata:
    """Load the wes_tumor_germline v1.2.1 example from grz-pydantic-models' own test resources."""
    with TEST_METADATA_PATH.open() as fh:
        return GrzSubmissionMetadata(**json.load(fh))


@pytest.fixture(
    params=[
        "sqlite",
        pytest.param(
            "postgresql",
            marks=pytest.mark.skipif(condition=which("pg_config") is None, reason="postgresql not detected"),
        ),
    ],
)
def db_test_connection(request: pytest.FixtureRequest):
    if request.param == "sqlite":
        tmpdir_factory: pytest.TempdirFactory = request.getfixturevalue("tmpdir_factory")
        db_dir = tmpdir_factory.mktemp("db")
        db_file = db_dir / "test.db"
        yield f"sqlite:///{str(db_file)}"
    elif request.param == "postgresql":
        postgresql: psycopg.Connection = request.getfixturevalue("postgresql")
        yield f"postgresql+psycopg://{postgresql.info.user}:@{postgresql.info.host}:{postgresql.info.port}/{postgresql.info.dbname}"
