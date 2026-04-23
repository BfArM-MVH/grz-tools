"""grz-check: Fast validation of sequencing files (FASTQ, BAM).

Native Python bindings to the Rust grz-check validator.

Examples:
    Recommended for production (large files, parallel validation) — mmap
    the file. Works for .fastq.gz too; the Rust side decompresses on the fly:

        >>> import mmap, grz_check
        >>> with open("huge.fastq.gz", "rb") as f, \\
        ...      mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as mm:
        ...     report = grz_check.validate_fastq(mm, min_mean_read_length=30)

    Simpler alternative — pass a path. Slower but fine for small files:

        >>> report = grz_check.validate_fastq("sample.fastq.gz", min_mean_read_length=30)

Any binary file-like object with ``.read()`` also works but is slower.
Writable buffers (``bytearray``, ``mmap`` opened for writing) are rejected.
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
        source: Path (str/PathLike), read-only buffer (bytes, memoryview, mmap
            ACCESS_READ — zero-copy), or a binary file-like object. For large
            files and parallel workloads, prefer ``mmap.mmap(f.fileno(), 0,
            access=mmap.ACCESS_READ)`` — works for .gz too; see module docs.
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
        r1: Read 1 as path (str/PathLike), read-only buffer, or file-like.
        r2: Read 2 as path (str/PathLike), read-only buffer, or file-like.
            For large files prefer mmap for both — see module docs.
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
