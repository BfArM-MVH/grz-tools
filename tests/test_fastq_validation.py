"""Tests for the fastq_validation module using ValidateOperation."""

from pathlib import Path

from grz_common.pipeline.operations import ValidateOperation


# Test case 1: Single end, line count not multiple of 4
def test_single_end_line_count_not_multiple_of_4():
    op = ValidateOperation()
    valid, errors, _stats = op.validate_fastq_file(
        Path("tests/mock_files/fastq_files_1000/single_end_failing.line_count.fastq.gz"), mean_read_length_threshold=75
    )

    assert not valid
    assert len(errors) >= 1
    assert "Number of lines is not a multiple of 4" in errors[0]


# Test case 2: Paired end, differing line numbers
def test_paired_end_differing_line_numbers():
    # Note: The specific logic to compare two files and raise a "Paired-end files have different read counts"
    # error has moved to the Submission worker.
    # Here we verify that the validator correctly reports line counts to enable that check.
    op = ValidateOperation()

    _valid1, _errors1, stats1 = op.validate_fastq_file(
        Path("tests/mock_files/fastq_files_1000/paired_end_failing_read1.differring_line_count.fastq.gz"),
        mean_read_length_threshold=75,
    )

    _valid2, _errors2, stats2 = op.validate_fastq_file(
        Path("tests/mock_files/fastq_files_1000/paired_end_failing_read2.differring_line_count.fastq.gz"),
        mean_read_length_threshold=75,
    )

    # Verify that we detected the mismatch in line counts
    assert stats1["line_count"] != stats2["line_count"]


# Test case 3: Paired end, all checks passed
def test_paired_end_all_checks_passed():
    op = ValidateOperation()

    valid1, errors1, stats1 = op.validate_fastq_file(
        Path("tests/mock_files/fastq_files_1000/paired_end_passing_read1.fastq.gz"),
        mean_read_length_threshold=75,
    )

    valid2, errors2, stats2 = op.validate_fastq_file(
        Path("tests/mock_files/fastq_files_1000/paired_end_passing_read2.fastq.gz"),
        mean_read_length_threshold=75,
    )

    assert valid1
    assert len(errors1) == 0
    assert valid2
    assert len(errors2) == 0

    # Verify line counts match
    assert stats1["line_count"] == stats2["line_count"]
