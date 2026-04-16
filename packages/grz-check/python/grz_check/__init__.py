"""grz-check: Fast validation of sequencing files (FASTQ, BAM)

This module provides native Python bindings to the Rust grz-check validator.

File-based validation (paths as strings):
    >>> import grz_check
    >>> report = grz_check.validate_fastq("sample.fastq.gz", min_mean_read_length=30)
    >>> print(f"Valid: {report.is_valid}, Records: {report.num_records}")

Stream-based validation (any file-like object with .read() method):
    >>> import gzip
    >>> import grz_check
    >>> with gzip.open("sample.fastq.gz", "rb") as f:
    ...     report = grz_check.validate_fastq(f, min_mean_read_length=30)
    >>> print(f"Valid: {report.is_valid}, Records: {report.num_records}")

Stream-based functions work with any object implementing the read trait:
    - file objects (open())
    - gzip streams (gzip.open())
    - BytesIO (io.BytesIO)
    - Custom stream classes with .read() method
"""

from __future__ import annotations

import os
import typing

from .grz_check import ValidationReport
from .grz_check import calculate_checksum as _calculate_checksum
from .grz_check import calculate_checksum_stream as _calculate_checksum_stream
from .grz_check import validate_bam as _validate_bam
from .grz_check import validate_bam_stream as _validate_bam_stream
from .grz_check import validate_fastq_paired_paths as _validate_fastq_paired_paths
from .grz_check import validate_fastq_paired_stream as _validate_fastq_paired_stream
from .grz_check import validate_fastq_single as _validate_fastq_single
from .grz_check import validate_fastq_single_stream as _validate_fastq_single_stream

__all__ = [
    "ValidationReport",
    "calculate_checksum",
    "validate_bam",
    "validate_fastq",
    "validate_fastq_paired",
]

__version__ = "0.2.1"


@typing.overload
def validate_fastq(source: str | os.PathLike, *, min_mean_read_length: int | None = None) -> ValidationReport: ...


@typing.overload
def validate_fastq(source: typing.BinaryIO, *, min_mean_read_length: int | None = None) -> ValidationReport: ...


def validate_fastq(source, *, min_mean_read_length=None):
    """Validate a single-end FASTQ file.

    Args:
        source: Path to the FASTQ file (str/PathLike) or a binary file-like object with .read() method.
        min_mean_read_length: Minimum mean read length required (>0), or None to skip check.

    Returns:
        ValidationReport with validation results.
    """
    if isinstance(source, (str, os.PathLike)):
        return _validate_fastq_single(str(source), min_mean_read_length=min_mean_read_length)
    return _validate_fastq_single_stream(source, min_mean_read_length=min_mean_read_length)


@typing.overload
def validate_fastq_paired(
    r1: str | os.PathLike, r2: str | os.PathLike, *, min_mean_read_length: int | None = None
) -> tuple[ValidationReport, ValidationReport]: ...


@typing.overload
def validate_fastq_paired(
    r1: typing.BinaryIO, r2: typing.BinaryIO, *, min_mean_read_length: int | None = None
) -> tuple[ValidationReport, ValidationReport]: ...


def validate_fastq_paired(r1, r2, *, min_mean_read_length=None):
    """Validate paired-end FASTQ files.

    Args:
        r1: Path (str/PathLike) or binary file-like object for read 1.
        r2: Path (str/PathLike) or binary file-like object for read 2.
        min_mean_read_length: Minimum mean read length required (>0), or None to skip check.

    Returns:
        Tuple of (r1_report, r2_report) ValidationReports.
    """
    if isinstance(r1, (str, os.PathLike)):
        return _validate_fastq_paired_paths(str(r1), str(r2), min_mean_read_length=min_mean_read_length)
    return _validate_fastq_paired_stream(r1, r2, min_mean_read_length=min_mean_read_length)


@typing.overload
def validate_bam(source: str | os.PathLike) -> ValidationReport: ...


@typing.overload
def validate_bam(source: typing.BinaryIO) -> ValidationReport: ...


def validate_bam(source):
    """Validate a BAM file.

    Args:
        source: Path to the BAM file (str/PathLike) or a binary file-like object with .read() method.

    Returns:
        ValidationReport with validation results.
    """
    if isinstance(source, (str, os.PathLike)):
        return _validate_bam(str(source))
    return _validate_bam_stream(source)


@typing.overload
def calculate_checksum(source: str | os.PathLike) -> str: ...


@typing.overload
def calculate_checksum(source: typing.BinaryIO) -> str: ...


def calculate_checksum(source):
    """Calculate SHA256 checksum.

    Args:
        source: Path to the file (str/PathLike) or a binary file-like object with .read() method.

    Returns:
        Hex-encoded SHA256 checksum string.
    """
    if isinstance(source, (str, os.PathLike)):
        return _calculate_checksum(str(source))
    return _calculate_checksum_stream(source)
