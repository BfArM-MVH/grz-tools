"""grz-check: Fast validation of sequencing files (FASTQ, BAM)

This module provides native Python bindings to the Rust grz-check validator.

File-based validation (paths as strings):
    >>> import grz_check
    >>> report = grz_check.validate_fastq_single("sample.fastq.gz", min_mean_read_length=30)
    >>> print(f"Valid: {report.is_valid}, Records: {report.num_records}")

Stream-based validation (any file-like object with .read() method):
    >>> import gzip
    >>> import grz_check
    >>> with gzip.open("sample.fastq.gz", "rb") as f:
    ...     report = grz_check.validate_fastq_single_stream(f, min_mean_read_length=30)
    >>> print(f"Valid: {report.is_valid}, Records: {report.num_records}")

Stream-based functions work with any object implementing the read trait:
    - file objects (open())
    - gzip streams (gzip.open())
    - BytesIO (io.BytesIO)
    - Custom stream classes with .read() method
"""

from .grz_check import (
    ValidationReport,
    calculate_checksum,
    calculate_checksum_stream,
    validate_bam,
    validate_bam_stream,
    validate_fastq_paired_paths,
    validate_fastq_paired_stream,
    validate_fastq_single,
    validate_fastq_single_stream,
)

__all__ = [
    "ValidationReport",
    "calculate_checksum",
    "calculate_checksum_stream",
    "validate_bam",
    "validate_bam_stream",
    "validate_fastq_paired_paths",
    "validate_fastq_paired_stream",
    "validate_fastq_single",
    "validate_fastq_single_stream",
]

__version__ = "0.2.1"
