import os

import pytest


@pytest.fixture(scope="session", autouse=True)
def mock_home(tmp_path_factory):
    # Create a temporary directory for the session
    temp_home = tmp_path_factory.mktemp("fake_home")

    # Set the environment variable
    os.environ["HOME"] = str(temp_home)
    os.environ["USERPROFILE"] = str(temp_home)  # For Windows compatibility

    return temp_home
