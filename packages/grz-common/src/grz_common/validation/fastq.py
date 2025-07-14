"""Validation of FASTQ files. Boils down to basic sanity checks such as line count and read length."""

import gzip
import logging
import os
import typing
from collections.abc import Generator
from contextlib import contextmanager
from gzip import GzipFile
from io import RawIOBase
from os import PathLike
from typing import TextIO

from tqdm.auto import tqdm

from ..constants import TQDM_SMOOTHING
from ..utils.io import TqdmIOWrapper

log = logging.getLogger(__name__)


def is_gzipped(file_path: str | PathLike) -> bool:
    """
    Check if a file is gzipped based on its extension.

    :param file_path: Path to the file
    :return: True if the file is gzipped, False otherwise
    """
    return str(file_path).endswith(".gz")


@contextmanager
def open_fastq(file_path: str | PathLike, progress=True) -> Generator[TextIO, None, None]:
    """
    Open a FASTQ file, handling both regular and gzipped formats.

    :param file_path: Path to the FASTQ file
    :param progress: Whether to show a progress bar
    :return: A file object opened in the appropriate mode (gzipped or plain text)
    """
    handle: TqdmIOWrapper | GzipFile | typing.BinaryIO | None = None
    with open(file_path, "rb") as fd:
        # Open FASTQ
        total_size = os.stat(file_path).st_size
        if progress:
            handle = TqdmIOWrapper(
                typing.cast(RawIOBase, fd),
                tqdm(
                    total=total_size,
                    desc="Reading FASTQ: ",
                    unit="B",
                    unit_scale=True,
                    # unit_divisor=1024,  # make use of standard units e.g. KB, MB, etc.
                    miniters=1,
                    smoothing=TQDM_SMOOTHING,
                ),
            )
        else:
            handle = fd

        if is_gzipped(file_path):
            # decompress
            with gzip.open(typing.cast(RawIOBase, handle), "rb") as decompressed_fd:
                yield typing.cast(TextIO, decompressed_fd)
        else:
            yield typing.cast(TextIO, handle)


def calculate_fastq_stats(file_path, expected_read_length: int | None = None) -> tuple[int, int | None]:
    """
    Calculate line number and read lengths in FASTQ file.
    :param file_path: Path to the FASTQ file
    :param expected_read_length: Expected read length (None if not known)
    :return: tuple with the following values:
      - Number of lines in the file
      - Observed read length
    """
    with open_fastq(file_path) as f:
        read_length_warned = False
        for line_number, line in enumerate(f):
            if (line_number % 4) == 1:
                # Sequence lines are every 4th line starting from the 2nd
                # Check if the read length is consistent
                read_length = len(line.strip())
                if expected_read_length is None:
                    expected_read_length = read_length
                elif (not read_length_warned) and (read_length != expected_read_length):
                    # For the time being, read length mismatch is downgraded to warning
                    log.warning(
                        f"Read length mismatch at line {line_number + 1}: "
                        f"expected {expected_read_length}, found {read_length}. "
                        "This will be an error in the future."
                        "Further mismatches in the same file won't be reported."
                    )
                    read_length_warned = True

    return (
        line_number + 1,  # enumerate starts indexing at 0
        expected_read_length,
    )


def validate_fastq_file(
    fastq_file: str | PathLike, expected_read_length: int | None = None
) -> tuple[int, int, list[str]] | tuple[int, int | None, list[str]]:
    """
    Validates a FASTQ file.

    :param fastq_file: Path to the fastq file
    :param expected_read_length: Expected read length (None if not known)
    :return: Tuple with the following fields:
      - number of lines
      - set of observed read lengths
      - list of errors found
    """
    errors = []
    try:
        if expected_read_length is not None:
            log.debug("%s: expecting read length of %s", fastq_file, expected_read_length)

        # Calculate the number of lines and read lengths
        num_lines, read_length = calculate_fastq_stats(fastq_file, expected_read_length=expected_read_length)
    except ValueError as e:
        if "Read length mismatch" in str(e):
            log.warning(f"{fastq_file}: {e}")
        else:
            errors.append(f"{fastq_file}: {e}")
        return -1, -1, errors

    # Check if the number of lines in a FASTQ file is a multiple of 4.
    if num_lines % 4 != 0:
        errors.append(f"{fastq_file}: Number of lines is not a multiple of 4! Found {num_lines} lines.")
    else:
        log.debug("%s: %s lines", fastq_file, num_lines)

    return num_lines, read_length, errors


def validate_single_end_reads(fastq_file: str | PathLike, expected_read_length: int | None = None) -> Generator[str]:
    """
    Validate a single-end FASTQ file.

    :param fastq_file: Path to the FASTQ file
    :param expected_read_length: Expected read length (None if not known)
    :return: Generator of errors, if any.
    """
    num_lines, read_lengths, errors = validate_fastq_file(fastq_file, expected_read_length=expected_read_length)
    yield from errors


def validate_paired_end_reads(
    fastq_file1: str | PathLike,
    fastq_file2: str | PathLike,
    expected_read_length: int | None = None,
) -> Generator[str]:
    """
    Validate two paired-end FASTQ files.

    :param fastq_file1: Path to the first FASTQ file (Read 1)
    :param fastq_file2: Path to the second FASTQ file (Read 2)
    :param expected_read_length: Expected read length (None if not known)
    :return: Generator of errors, if any.
    """
    num_lines_file1, read_lengths_file1, errors_file1 = validate_fastq_file(
        fastq_file1, expected_read_length=expected_read_length
    )
    yield from errors_file1
    num_lines_file2, read_lengths_file2, errors_file2 = validate_fastq_file(
        fastq_file2, expected_read_length=expected_read_length
    )
    yield from errors_file2

    if num_lines_file1 != num_lines_file2:
        yield f"Paired-end files have different read counts: '{fastq_file1}' ({num_lines_file1}) and '{fastq_file2}' ({num_lines_file2})!"
