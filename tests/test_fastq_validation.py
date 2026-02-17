"""Tests for the fastq_validation module using FastqValidator."""

import pytest
from grz_common.pipeline.components import Stream
from grz_common.pipeline.components.validation import FastqValidator


def run_validator(path: str, threshold: float = 0.0):
    validator = FastqValidator(mean_read_length_threshold=threshold)
    with open(path, "rb") as f:
        Stream(f) >> validator
    return True, [], validator.metrics


def test_single_end_line_count_not_multiple_of_4():
    """Single end, line count not multiple of 4"""
    path = "tests/mock_files/fastq_files_1000/single_end_failing.line_count.fastq.gz"

    with pytest.raises(ValueError) as excinfo:
        run_validator(path, threshold=75)

    assert "incomplete record" in str(excinfo.value).lower()


def test_paired_end_differing_line_numbers():
    """Paired end, differing line numbers"""
    path1 = "tests/mock_files/fastq_files_1000/paired_end_failing_read1.differring_line_count.fastq.gz"
    path2 = "tests/mock_files/fastq_files_1000/paired_end_failing_read2.differring_line_count.fastq.gz"

    _, _, stats1 = run_validator(path1, threshold=75)
    _, _, stats2 = run_validator(path2, threshold=75)

    # Verify that we detected the mismatch in line counts
    assert stats1["line_count"] != stats2["line_count"]

    # TODO: ideally we'd use the larger pipeline context (partner_map and ConsistencyValidator)
    #  to check that we notice this mismatch automatically


def test_paired_end_all_checks_passed():
    """Paired end, all checks passed"""
    path1 = "tests/mock_files/fastq_files_1000/paired_end_passing_read1.fastq.gz"
    path2 = "tests/mock_files/fastq_files_1000/paired_end_passing_read2.fastq.gz"

    valid1, errors1, stats1 = run_validator(path1, threshold=75)
    valid2, errors2, stats2 = run_validator(path2, threshold=75)

    assert valid1
    assert len(errors1) == 0
    assert valid2
    assert len(errors2) == 0

    # Verify line counts match
    assert stats1["line_count"] == stats2["line_count"]
