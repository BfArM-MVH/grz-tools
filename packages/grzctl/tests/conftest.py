import os

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

# Rich decides whether to emit ANSI colour when a Console is constructed, which for the CLI modules
# happens at import. TTY_COMPATIBLE=0 outranks both FORCE_COLOR and isatty detection, so captured
# output stays plain whatever the developer's shell exports. Set at import: a fixture would run
# after the test modules have already built their consoles.
os.environ["TTY_COMPATIBLE"] = "0"


@pytest.fixture(scope="session", autouse=True)
def mock_home(tmp_path_factory):
    # Create a temporary directory for the session
    temp_home = tmp_path_factory.mktemp("fake_home")

    # Set the environment variable
    os.environ["HOME"] = str(temp_home)
    os.environ["USERPROFILE"] = str(temp_home)  # For Windows compatibility

    return temp_home
