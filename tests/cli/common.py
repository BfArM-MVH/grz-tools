"""
Common methods/fixtures for CLI tests
"""

import shutil
from pathlib import Path

import pytest
import requests
from urllib3.exceptions import MaxRetryError, NewConnectionError

SUBMISSION_DIR = Path("tests/mock_files/submissions/valid_submission")


def copy_submission(working_dir_path: Path, *parts: str, source: Path = SUBMISSION_DIR) -> None:
    """Copy the named parts of the example submission into the working directory.

    With no parts given, the entire source directory is copied instead.
    """
    if not parts:
        shutil.copytree(source, working_dir_path, dirs_exist_ok=True)
        return
    for part in parts:
        shutil.copytree(source / part, working_dir_path / part, dirs_exist_ok=True)


def refused_connection(url: str) -> requests.ConnectionError:
    """What requests raises when the server of ``url`` refuses the connection, so that nothing was sent."""
    return requests.ConnectionError(MaxRetryError(None, url, NewConnectionError(None, "Connection refused")))


@pytest.fixture
def working_dir(tmpdir_factory: pytest.TempdirFactory):
    """Create temporary folder for the session"""
    datadir = tmpdir_factory.mktemp("submission")
    return datadir


@pytest.fixture
def working_dir_path(working_dir) -> Path:
    return Path(working_dir.strpath)
