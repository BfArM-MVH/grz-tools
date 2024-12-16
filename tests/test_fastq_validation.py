"""Tests for the fastq_validation module."""

from grz_cli.fastq_validation import (
    validate_paired_end_reads,
    validate_single_end_reads,
)


# Test case 1: Single end, line count not multiple of 4
def test_single_end_line_count_not_multiple_of_4():
    errors = list(validate_single_end_reads("tests/mock_files/fastq_files_1000/single_end_failing.line_count.fastq.gz"))
    assert len(errors) == 1
    assert "Number of lines is not a multiple of 4" in errors[0]


# Test case 2: Single end, inconsistent read length
def test_single_end_inconsistent_read_length():
    errors = list(
        validate_single_end_reads(
            "tests/mock_files/fastq_files_1000/single_end_failing.inconsistent_read_length.fastq.gz"
        )
    )
    assert len(errors) == 1
    assert "Inconsistent read lengths" in errors[0]


# Test case 3: Paired end, differing line numbers
def test_paired_end_differing_line_numbers():
    errors = list(
        validate_paired_end_reads(
            "tests/mock_files/fastq_files_1000/paired_end_failing_read1.differring_line_count.fastq.gz",
            "tests/mock_files/fastq_files_1000/paired_end_failing_read2.differring_line_count.fastq.gz",
        )
    )
    assert len(errors) == 1
    assert "Paired-end files have different read counts" in errors[0]


# Test case 4: Paired end, all checks passed
def test_paired_end_all_checks_passed():
    errors = list(
        validate_paired_end_reads(
            "tests/mock_files/fastq_files_1000/paired_end_passing_read1.fastq.gz",
            "tests/mock_files/fastq_files_1000/paired_end_passing_read2.fastq.gz",
        )
    )
    assert len(errors) == 0
