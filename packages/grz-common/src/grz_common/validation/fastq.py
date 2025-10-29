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
from pathlib import Path
from typing import TextIO

from tqdm.auto import tqdm

from ..constants import TQDM_DEFAULTS
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
    file_name = Path(file_path).name
    with open(file_path, "rb") as fd:
        # Open FASTQ
        total_size = os.stat(file_path).st_size
        if progress:
            handle = TqdmIOWrapper(
                typing.cast(RawIOBase, fd),
                tqdm(total=total_size, desc="FASTQ   ", postfix=f"{file_name}", **TQDM_DEFAULTS),  # type: ignore[call-overload]
            )
        else:
            handle = fd

        if is_gzipped(file_path):
            # decompress
            with gzip.open(typing.cast(RawIOBase, handle), "rb") as decompressed_fd:
                yield typing.cast(TextIO, decompressed_fd)
        else:
            yield typing.cast(TextIO, handle)


def calculate_fastq_stats(file_path) -> tuple[int, float]:
    """
    Calculate line number and read lengths in FASTQ file.

    :param file_path: Path to the FASTQ file
    :return: tuple with the following values:
      - Number of lines in the file
      - Observed mean read length
    """
    total_read_length = 0
    total_reads = 0
    with open_fastq(file_path) as f:
        for line_number, line in enumerate(f):
            if (line_number % 4) == 1:
                # Sequence lines are every 4th line starting from the 2nd
                total_read_length += len(line.strip())
                total_reads += 1

    return (
        line_number + 1,  # enumerate starts indexing at 0
        total_read_length / total_reads,
    )


def validate_fastq_file(fastq_file: str | PathLike, mean_read_length_threshold: int) -> tuple[int, list[str]]:
    """
    Validates a FASTQ file.

    :param fastq_file: Path to the fastq file
    :param mean_read_length_threshold: Exclusive minimum mean read length
    :return: Tuple with the following fields:
      - number of lines
      - set of observed read lengths
      - list of errors found
    """
    errors = []
    try:
        # Calculate the number of lines and read lengths
        num_lines, mean_read_length = calculate_fastq_stats(fastq_file)
    except ValueError as e:
        errors.append(f"{fastq_file}: {e}")
        return -1, errors

    # Check if the number of lines in a FASTQ file is a multiple of 4.
    if num_lines % 4 != 0:
        errors.append(f"{fastq_file}: Number of lines is not a multiple of 4! Found {num_lines} lines.")
    else:
        log.debug("%s: %s lines", fastq_file, num_lines)

    if mean_read_length <= mean_read_length_threshold:
        raise ValueError(
            f"Mean read length must be > {mean_read_length_threshold}bp, calculated {mean_read_length:.2f}."
        )

    return num_lines, errors


def validate_single_end_reads(fastq_file: str | PathLike, mean_read_length_threshold: int) -> Generator[str]:
    """
    Validate a single-end FASTQ file.

    :param fastq_file: Path to the FASTQ file
    :return: Generator of errors, if any.
    :param mean_read_length_threshold: Exclusive minimum mean read length
    """
    num_lines, errors = validate_fastq_file(fastq_file, mean_read_length_threshold=mean_read_length_threshold)
    yield from errors


def validate_paired_end_reads(
    fastq_file1: str | PathLike, fastq_file2: str | PathLike, mean_read_length_threshold: int
) -> Generator[str]:
    """
    Validate two paired-end FASTQ files.

    :param fastq_file1: Path to the first FASTQ file (Read 1)
    :param fastq_file2: Path to the second FASTQ file (Read 2)
    :param mean_read_length_threshold: Exclusive minimum mean read length
    :return: Generator of errors, if any.
    """
    num_lines_file1, errors_file1 = validate_fastq_file(
        fastq_file1, mean_read_length_threshold=mean_read_length_threshold
    )
    yield from errors_file1
    num_lines_file2, errors_file2 = validate_fastq_file(
        fastq_file2, mean_read_length_threshold=mean_read_length_threshold
    )
    yield from errors_file2

    if num_lines_file1 != num_lines_file2:
        yield f"Paired-end files have different read counts: '{fastq_file1}' ({num_lines_file1}) and '{fastq_file2}' ({num_lines_file2})!"
