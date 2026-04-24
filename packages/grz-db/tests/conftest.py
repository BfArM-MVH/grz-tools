"""Shared fixtures for grz-db tests."""

import importlib.util
import json
from pathlib import Path
from shutil import which

import psycopg
import pytest
from grz_pydantic_models.submission.metadata import GrzSubmissionMetadata

# Resolve the grz-pydantic-models package root, then navigate to its test resources.
# spec.origin  →  .../packages/grz-pydantic-models/src/grz_pydantic_models/__init__.py
# parents[2]   →  .../packages/grz-pydantic-models/
_GRZ_PYDANTIC_MODELS_SPEC = importlib.util.find_spec("grz_pydantic_models")
assert _GRZ_PYDANTIC_MODELS_SPEC is not None and _GRZ_PYDANTIC_MODELS_SPEC.origin is not None
_GRZ_PYDANTIC_MODELS_ROOT = Path(_GRZ_PYDANTIC_MODELS_SPEC.origin).parents[2]
_MOCK_METADATA_PATH = _GRZ_PYDANTIC_MODELS_ROOT / "tests/resources/example_metadata/wes_tumor_germline/v1.2.1.json"


@pytest.fixture(scope="session")
def metadata() -> GrzSubmissionMetadata:
    """Load the wes_tumor_germline v1.2.1 example from grz-pydantic-models' own test resources."""
    with _MOCK_METADATA_PATH.open() as fh:
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
