"""Tests for the fastq_validation module."""

import grz_check


# Test case 1: Single end, line count not multiple of 4
def test_single_end_line_count_not_multiple_of_4():
    report = grz_check.validate_fastq(
        "tests/mock_files/fastq_files_1000/single_end_failing.line_count.fastq.gz", min_mean_read_length=75
    )

    assert not report.is_valid
    assert len(report.errors) >= 1
    assert any(
        "failed to fill whole buffer" in e.lower() or "failed to parse record" in e.lower() for e in report.errors
    )


# Test case 2: Paired end, differing line numbers
def test_paired_end_differing_line_numbers():
    r1_report, r2_report = grz_check.validate_fastq_paired(
        "tests/mock_files/fastq_files_1000/paired_end_failing_read1.differring_line_count.fastq.gz",
        "tests/mock_files/fastq_files_1000/paired_end_failing_read2.differring_line_count.fastq.gz",
        min_mean_read_length=75,
    )
    print(r1_report.errors, r1_report.warnings, r1_report.num_records, r1_report.is_valid, r1_report.mean_read_length)
    print(r2_report.errors, r2_report.warnings, r2_report.num_records, r2_report.is_valid, r2_report.mean_read_length)

    errors = list(r1_report.errors) + list(r2_report.errors)

    # At least one of the reports should be marked invalid
    assert not r1_report.is_valid or not r2_report.is_valid
    assert len(errors) >= 1
    assert any("Mismatched read counts" in e for e in errors)


# Test case 3: Paired end, all checks passed
def test_paired_end_all_checks_passed():
    r1_report, r2_report = grz_check.validate_fastq_paired(
        "tests/mock_files/fastq_files_1000/paired_end_passing_read1.fastq.gz",
        "tests/mock_files/fastq_files_1000/paired_end_passing_read2.fastq.gz",
        min_mean_read_length=75,
    )

    assert r1_report.is_valid
    assert len(r1_report.errors) == 0

    assert r2_report.is_valid
    assert len(r2_report.errors) == 0
