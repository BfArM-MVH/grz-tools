"""Pytest fixtures for provisioning :class:`~grz_db.models.submission.SubmissionDb` databases.

Requires the ``testing`` extra. Import the fixtures into a ``conftest.py``::

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

The underscore-prefixed session fixtures are resolved by name at runtime, so they have to be
imported too. Importing rather than declaring ``pytest_plugins`` is deliberate twice over: pytest
only accepts that declaration in a conftest at the rootdir, and a rootdir declaration would load
this module for every package in the workspace, including those whose environments have no
pytest-postgresql. No ``pytest11`` entry point is registered for the same reason.

Each database is cloned from a template migrated once per session, rather than migrated per test.
The template is a real Alembic run because the partial unique indexes live in migrations rather
than in the models, so a schema built from ``SQLModel.metadata.create_all`` would silently lack
them.
"""

import itertools
import shutil
from collections.abc import Iterator
from pathlib import Path

import psycopg
import pytest
import sqlalchemy
from cryptography.hazmat.primitives import serialization
from cryptography.hazmat.primitives.asymmetric import ed25519
from pytest_postgresql import factories

from grz_db.models.author import Author
from grz_db.models.submission import SubmissionDb

__all__ = [
    "db",
    "db_test_connection",
    "migrated_db_connection",
    "postgresql_proc",
    "test_author",
]

_PG_TEMPLATE_DBNAME = "grz_db_migrated_template"

_BACKEND_PARAMS = [
    "sqlite",
    pytest.param(
        "postgresql",
        marks=pytest.mark.skipif(condition=shutil.which("pg_config") is None, reason="postgresql not detected"),
    ),
]

# The cluster lives and dies with the test session, so crash safety protects nothing while an
# fsync per commit costs most of the run.
postgresql_proc = factories.postgresql_proc(
    postgres_options="-c fsync=off -c synchronous_commit=off -c full_page_writes=off",
)


def _migrate(db_url: str) -> None:
    """Bring *db_url* to the latest schema and release its connection pool.

    :param db_url: SQLAlchemy URL of an empty database.
    """
    submission_db = SubmissionDb(db_url=db_url, author=None)
    try:
        submission_db.initialize_schema()
    finally:
        submission_db.engine.dispose()


def _pg_url(proc, dbname: str) -> str:
    """Build a SQLAlchemy URL for *dbname* on the cluster *proc* runs."""
    return f"postgresql+psycopg://{proc.user}:@{proc.host}:{proc.port}/{dbname}"


def _pg_admin_connect(proc) -> psycopg.Connection:
    """Connect to the maintenance database, from which others can be created and dropped."""
    return psycopg.connect(host=proc.host, port=proc.port, user=proc.user, dbname="postgres", autocommit=True)


def skip_sqlite_fsync(engine: sqlalchemy.Engine) -> None:
    """Stop SQLite flushing to disk on every commit.

    These databases are discarded when the test ends, so the durability that ``synchronous=FULL``
    buys protects nothing and costs an fsync per transaction. Tests that commit in a loop spend
    almost all of their time there.

    :param engine: Engine to configure. Engines on other backends are left alone.
    """
    if engine.dialect.name != "sqlite":
        return

    @sqlalchemy.event.listens_for(engine, "connect")
    def _set_pragma(dbapi_connection, _connection_record):
        dbapi_connection.execute("PRAGMA synchronous=OFF")


@pytest.fixture(params=_BACKEND_PARAMS)
def db_backend(request: pytest.FixtureRequest) -> str:
    """The backend under test, one of ``sqlite`` or ``postgresql``.

    Fixtures that provision a database take their parametrization from this one, so a test
    requesting several of them still runs once per backend rather than once per combination.
    """
    return str(request.param)


@pytest.fixture
def db_test_connection(
    db_backend: str, request: pytest.FixtureRequest, tmp_path_factory: pytest.TempPathFactory
) -> Iterator[str]:
    """Database URL of an empty database, one per supported backend.

    For tests that drive the migrations themselves. Anything wanting the current schema should use
    :func:`migrated_db_connection`, which is far cheaper per test.
    """
    if db_backend == "sqlite":
        yield f"sqlite:///{tmp_path_factory.mktemp('db') / 'test.db'}"
    else:
        postgresql: psycopg.Connection = request.getfixturevalue("postgresql")
        info = postgresql.info
        yield f"postgresql+psycopg://{info.user}:@{info.host}:{info.port}/{info.dbname}"


@pytest.fixture(scope="session")
def _migrated_sqlite_template(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """A migrated SQLite file that each test copies."""
    template = tmp_path_factory.mktemp("grz-db-template") / "template.db"
    _migrate(f"sqlite:///{template}")
    return template


@pytest.fixture(scope="session")
def _migrated_postgresql_template(request: pytest.FixtureRequest) -> str:
    """Name of a migrated database that each test clones."""
    proc = request.getfixturevalue("postgresql_proc")
    with _pg_admin_connect(proc) as conn:
        conn.execute(f'DROP DATABASE IF EXISTS "{_PG_TEMPLATE_DBNAME}"')
        conn.execute(f'CREATE DATABASE "{_PG_TEMPLATE_DBNAME}"')
    _migrate(_pg_url(proc, _PG_TEMPLATE_DBNAME))
    return _PG_TEMPLATE_DBNAME


@pytest.fixture(scope="session")
def _pg_test_dbnames() -> Iterator[int]:
    """Source of database names unique within the session."""
    return itertools.count()


@pytest.fixture
def migrated_db_connection(
    db_backend: str, request: pytest.FixtureRequest, tmp_path_factory: pytest.TempPathFactory
) -> Iterator[str]:
    """Database URL of a database already on the latest schema, one per supported backend.

    The schema is built once per session and cloned per test, so a test pays a file copy or a
    ``CREATE DATABASE ... TEMPLATE`` rather than a full Alembic run.
    """
    if db_backend == "sqlite":
        template: Path = request.getfixturevalue("_migrated_sqlite_template")
        db_file = tmp_path_factory.mktemp("db") / "test.db"
        shutil.copyfile(template, db_file)
        yield f"sqlite:///{db_file}"
    else:
        template_dbname: str = request.getfixturevalue("_migrated_postgresql_template")
        proc = request.getfixturevalue("postgresql_proc")
        dbname = f"grz_db_test_{next(request.getfixturevalue('_pg_test_dbnames'))}"
        with _pg_admin_connect(proc) as conn:
            conn.execute(f'CREATE DATABASE "{dbname}" TEMPLATE "{template_dbname}"')
        try:
            yield _pg_url(proc, dbname)
        finally:
            # FORCE: a test that leaves a connection open would otherwise block the drop.
            with _pg_admin_connect(proc) as conn:
                conn.execute(f'DROP DATABASE IF EXISTS "{dbname}" WITH (FORCE)')


@pytest.fixture
def test_author() -> Author:
    """An author signing as ``alice``, with a throwaway key."""
    key = ed25519.Ed25519PrivateKey.generate()
    private_key_bytes = key.private_bytes(
        encoding=serialization.Encoding.PEM,
        format=serialization.PrivateFormat.OpenSSH,
        encryption_algorithm=serialization.NoEncryption(),
    )
    return Author(name="alice", private_key_bytes=private_key_bytes, private_key_passphrase="")


@pytest.fixture
def db(migrated_db_connection: str, test_author: Author) -> Iterator[SubmissionDb]:
    """A ``SubmissionDb`` on the latest schema, one per supported backend."""
    submission_db = SubmissionDb(db_url=migrated_db_connection, author=test_author, debug=False)
    skip_sqlite_fsync(submission_db.engine)
    try:
        yield submission_db
    finally:
        submission_db.engine.dispose()
