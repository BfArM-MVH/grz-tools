"""Tests for ``Worker.validate``'s error handling."""

import grz_common.exceptions as grzexc
import pytest
from grz_common.models.identifiers import IdentifiersModel
from grz_common.workers.submission import Submission
from grz_common.workers.worker import Worker


@pytest.fixture
def worker(submission_metadata_dir, tmp_path) -> Worker:
    return Worker(
        metadata_dir=submission_metadata_dir,
        files_dir=tmp_path / "files",
        log_dir=tmp_path / "logs",
        encrypted_files_dir=tmp_path / "encrypted_files",
    )


def test_a_failed_file_validation_is_not_wrapped_a_second_time(worker, monkeypatch):
    """``validate`` already raises ``SubmissionValidationError`` for a failed file validation.

    Its own ``except Exception`` must not wrap that error again, or the message and the log
    both double up.
    """
    monkeypatch.setattr(Submission, "validate_files", lambda self, **kwargs: iter(["a file failed"]))

    with pytest.raises(grzexc.SubmissionValidationError) as excinfo:
        worker.validate(identifiers=IdentifiersModel(grz="GRZK00007", le="260914050"))

    assert not str(excinfo.value).startswith("Validation failed due to an error")
