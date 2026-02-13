"""Tests for the Submission worker validation logic."""

from pathlib import Path
from unittest.mock import MagicMock

import pytest
from grz_common.pipeline.operations import ValidateOperation
from grz_common.workers.submission import Submission
from grz_pydantic_models.submission.metadata.v1 import ReadOrder


@pytest.fixture
def mock_submission():
    """Create a Submission instance with mocked init."""
    # bypass __init__ to avoid loading metadata file from disk
    submission = Submission.__new__(Submission)
    submission.files_dir = Path(".")

    return submission


@pytest.fixture
def mock_logger():
    """Mock FileProgressLogger."""
    logger = MagicMock()
    # Always return None for get_state to force re-validation
    logger.get_state.return_value = None
    return logger


def test_paired_end_fallback_differing_line_counts(mock_submission, mock_logger):
    """
    Test that paired-end validation detects differing line counts.
    """
    r1_path = "tests/mock_files/fastq_files_1000/paired_end_failing_read1.differring_line_count.fastq.gz"
    r2_path = "tests/mock_files/fastq_files_1000/paired_end_failing_read2.differring_line_count.fastq.gz"

    f1 = MagicMock()
    f1.file_path = r1_path
    f1.read_order = ReadOrder.r1
    f1.flowcell_id = "FC1"
    f1.lane_id = 1

    f2 = MagicMock()
    f2.file_path = r2_path
    f2.read_order = ReadOrder.r2
    f2.flowcell_id = "FC1"
    f2.lane_id = 1

    fastq_files = [f1, f2]

    thresholds = MagicMock()
    thresholds.mean_read_length = 0

    errors = list(mock_submission._validate_paired_end_fallback(fastq_files, mock_logger, thresholds))

    assert len(errors) >= 1
    assert any("Paired-end files have different read counts" in e for e in errors)


def test_paired_end_fallback_success(mock_submission, mock_logger):
    """Test that valid paired-end files pass the fallback validation logic."""
    r1_path = "tests/mock_files/fastq_files_1000/paired_end_passing_read1.fastq.gz"
    r2_path = "tests/mock_files/fastq_files_1000/paired_end_passing_read2.fastq.gz"

    f1 = MagicMock()
    f1.file_path = r1_path
    f1.read_order = ReadOrder.r1
    f1.flowcell_id = "FC1"
    f1.lane_id = 1

    f2 = MagicMock()
    f2.file_path = r2_path
    f2.read_order = ReadOrder.r2
    f2.flowcell_id = "FC1"
    f2.lane_id = 1

    fastq_files = [f1, f2]

    thresholds = MagicMock()
    thresholds.mean_read_length = 0

    errors = list(mock_submission._validate_paired_end_fallback(fastq_files, mock_logger, thresholds))

    assert len(errors) == 0
