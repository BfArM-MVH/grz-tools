"""grz-check: Fast validation of sequencing files (FASTQ, BAM)

This module provides native Python bindings to the Rust grz-check validator.

Example:
    >>> import grz_check
    >>> report = grz_check.validate_fastq_single("sample.fastq.gz", min_mean_read_length=30)
    >>> print(f"Valid: {report.is_valid}, Records: {report.num_records}")
"""

from .grz_check import (
    ValidationReport,
    validate_fastq_single,
    validate_fastq_paired,
    validate_bam,
    calculate_checksum,
)

__all__ = [
    "ValidationReport",
    "validate_fastq_single",
    "validate_fastq_paired",
    "validate_bam",
    "calculate_checksum",
]

__version__ = "0.2.1"
