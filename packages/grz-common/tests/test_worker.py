"""``Worker`` sets the mode of the directories it writes into only if it is given a mode."""

import os
import stat
from pathlib import Path

import pytest
from grz_common.workers.worker import Worker


def _mode(path: Path) -> int:
    return stat.S_IMODE(path.stat().st_mode)


def _worker(submission_dir: Path, directory_mode: int | None) -> Worker:
    return Worker(
        metadata_dir=submission_dir / "metadata",
        files_dir=submission_dir / "files",
        log_dir=submission_dir / "logs",
        encrypted_files_dir=submission_dir / "encrypted_files",
        directory_mode=directory_mode,
    )


@pytest.fixture
def submission_dir(tmp_path: Path) -> Path:
    submission_dir = tmp_path / "submission"
    submission_dir.mkdir()
    submission_dir.chmod(0o755)
    return submission_dir


def test_worker_without_directory_mode_creates_the_log_directory_under_the_umask(submission_dir: Path):
    """grz-cli passes the directories of a submitter, so the umask decides their mode, as before."""
    previous_umask = os.umask(0o022)
    try:
        _worker(submission_dir, directory_mode=None)
    finally:
        os.umask(previous_umask)

    assert _mode(submission_dir / "logs") == 0o750


def test_worker_without_directory_mode_keeps_the_mode_of_an_existing_log_directory(submission_dir: Path):
    (submission_dir / "logs").mkdir()
    (submission_dir / "logs").chmod(0o755)

    _worker(submission_dir, directory_mode=None)

    assert _mode(submission_dir / "logs") == 0o755


def test_worker_with_directory_mode_sets_it_on_an_existing_log_directory(submission_dir: Path):
    (submission_dir / "logs").mkdir()
    (submission_dir / "logs").chmod(0o755)

    _worker(submission_dir, directory_mode=0o700)

    assert _mode(submission_dir / "logs") == 0o700
    assert _mode(submission_dir) == 0o755
