from shutil import which

import psycopg
import pytest


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
