"""Tests for the file_progress_logging module."""

import hashlib
import json
from pathlib import Path

import pytest
from grz_pydantic_models.v1_1_1.metadata import File as SubmissionFileMetadata

from grz_cli.file_operations import read_multiple_json
from grz_cli.progress_logging import FileProgressLogger


@pytest.fixture(scope="function")
def test_dir(tmpdir_factory: pytest.TempdirFactory):
    datadir = tmpdir_factory.mktemp("test")
    return datadir


@pytest.fixture
def test_dir_path(test_dir) -> Path:
    return Path(test_dir.strpath)


@pytest.fixture
def temp_log_file_path(test_dir_path) -> Path:
    """Fixture to create a temporary log file."""
    temp_file = test_dir_path / "log.mjson"
    return temp_file


@pytest.fixture
def temp_data_file_path(test_dir_path) -> Path:
    """Fixture to create a temporary file to track."""
    temp_file = test_dir_path / "data.bed"
    with open(temp_file, "w") as fd:
        fd.write("asdf")

    return temp_file


@pytest.fixture
def temp_file_metadata_dict(temp_data_file_path):
    with open(temp_data_file_path, "rb") as f:
        sha256sum = hashlib.sha256(f.read()).hexdigest()

    stat = temp_data_file_path.stat()

    return {
        "filePath": str(temp_data_file_path),
        "fileType": "bed",
        "fileChecksum": sha256sum,
        "fileSizeInBytes": stat.st_size,
        "checksumType": "sha256",
    }


@pytest.fixture
def temp_data_file_metadata(temp_file_metadata_dict):
    return SubmissionFileMetadata(**temp_file_metadata_dict)


@pytest.fixture
def logger(temp_log_file_path):
    """Fixture to create a FileProgressLogger instance with sample fields."""
    return FileProgressLogger(temp_log_file_path)


def test_set_and_get_state(
    logger: FileProgressLogger,
    temp_data_file_path: Path,
    temp_file_metadata_dict: dict,
    temp_data_file_metadata: SubmissionFileMetadata,
):
    """Test that setting and getting the state works."""
    # Set state for a file
    state = {"progress": 50, "status": "in_progress"}
    logger.set_state(
        file_path=temp_data_file_path,
        file_metadata=SubmissionFileMetadata(**dict(temp_file_metadata_dict)),
        state=state,
    )

    # Retrieve the state
    retrieved_state = logger.get_state(temp_data_file_path, temp_data_file_metadata)
    assert retrieved_state == state


def test_get_state_file_not_tracked(
    logger: FileProgressLogger,
    temp_data_file_path: Path,
    temp_data_file_metadata: SubmissionFileMetadata,
):
    """Test that get_state returns None for untracked files."""
    state = {"progress": 75, "status": "completed"}
    logger.set_state(
        file_path=temp_data_file_path,
        file_metadata=temp_data_file_metadata,
        state=state,
    )
    # Retrieve the state
    retrieved_state = logger.get_state(
        temp_data_file_path,
        file_metadata=SubmissionFileMetadata(
            **{
                "filePath": "foo",
                "fileType": "bed",
                "fileChecksum": "b5bb9d8014a0f9b1d61e21e796d78dccdf1352f23cd32812f4850b878ae4944c",
                "fileSizeInBytes": 0,
            }
        ),
    )
    assert retrieved_state is None


def test_persist_state_to_json(
    logger,
    temp_log_file_path,
    temp_data_file_path,
    temp_data_file_metadata,
):
    """Test that file state is persisted to json."""
    state = {"progress": 75, "status": "completed"}
    logger.set_state(
        file_path=temp_data_file_path,
        file_metadata=temp_data_file_metadata,
        state=state,
    )

    # Check that the state was appended to the json file
    with open(temp_log_file_path) as f:
        entries = list(read_multiple_json(f))
        assert len(entries) == 1
        assert entries[0] == {
            "file_path": str(temp_data_file_path),
            "modification_time": temp_data_file_path.stat().st_mtime,
            "size": temp_data_file_path.stat().st_size,
            "metadata": temp_data_file_metadata.model_dump(by_alias=True),
            "state": {"progress": 75, "status": "completed"},
        }


def test_get_index(mocker, logger, temp_data_file_path):
    """Test that file index is generated correctly."""
    # Mock the stat result to control modification time
    real_st_mode = temp_data_file_path.stat().st_mode

    mock_stat = mocker.patch("pathlib.Path.stat")
    mock_stat.return_value.st_mtime = 1234567890.0
    mock_stat.return_value.st_size = 1234
    mock_stat.return_value.st_mode = real_st_mode

    index = logger._get_index(temp_data_file_path)
    assert index == (str(temp_data_file_path), 1234567890.0, 1234)


def test_read_existing_log(temp_log_file_path, temp_data_file_path, temp_data_file_metadata):
    """Test that file states are correctly read from an existing json log."""
    # Manually write a row to the log file
    with open(temp_log_file_path, "w", newline="") as fd:
        json.dump(
            {
                "file_path": str(temp_data_file_path),
                "modification_time": temp_data_file_path.stat().st_mtime,
                "size": temp_data_file_path.stat().st_size,
                "metadata": temp_data_file_metadata.model_dump(by_alias=True),
                "state": {"progress": 0, "status": "started"},
            },
            fd,
        )
        fd.write("\n")
        json.dump(
            {
                "file_path": str(temp_data_file_path),
                "modification_time": temp_data_file_path.stat().st_mtime,
                "size": temp_data_file_path.stat().st_size,
                "metadata": temp_data_file_metadata.model_dump(by_alias=True),
                "state": {"progress": 25, "status": "in-progress"},
            },
            fd,
        )
        fd.write("\n")

    # Reload the logger to test reading from the log
    logger = FileProgressLogger(temp_log_file_path)

    # cleanup
    logger.cleanup(keep=[(temp_data_file_path, temp_data_file_metadata)])

    # Verify the state is read correctly
    retrieved_state = logger.get_state(temp_data_file_path, temp_data_file_metadata)
    assert retrieved_state is not None
    assert retrieved_state == {"progress": 25, "status": "in-progress"}
