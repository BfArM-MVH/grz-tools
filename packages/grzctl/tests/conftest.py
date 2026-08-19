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

# Rich decides whether to emit ANSI colour when a Console is constructed, and for the CLI modules
# that happens at import time. TTY_COMPATIBLE=0 outranks both FORCE_COLOR and isatty detection, so
# captured output stays plain no matter what the developer's shell exports. This is set here, at
# import time, because a fixture would only run after the test modules have already built their
# consoles.
os.environ["TTY_COMPATIBLE"] = "0"


@pytest.fixture(scope="session", autouse=True)
def mock_home(tmp_path_factory):
    # Create a temporary directory for the session
    temp_home = tmp_path_factory.mktemp("fake_home")

    # Set the environment variable
    os.environ["HOME"] = str(temp_home)
    os.environ["USERPROFILE"] = str(temp_home)  # For Windows compatibility

    return temp_home
