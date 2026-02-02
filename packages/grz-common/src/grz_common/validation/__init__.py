"""Validation package for grz-cli."""

import logging
import subprocess
from pathlib import Path


class UserInterruptException(Exception):
    """Raised when an interrupt is triggered via ctrl+c."""

    pass


# Gzip magic bytes
GZIP_MAGIC = b"\x1f\x8b"


def check_gzip_magic_bytes(file_path: Path) -> str | None:
    """
    Check if a file with .gz extension has valid gzip magic bytes.

    :param file_path: Path to the file to check
    :return: Error message if validation fails, None otherwise
    """
    # Only check files with .gz extension
    if not str(file_path).endswith(".gz"):
        return None

    try:
        with open(file_path, "rb") as f:
            magic_bytes = f.read(2)

        if len(magic_bytes) < 2:
            return f"File has .gz extension but is too small to contain gzip header (size: {len(magic_bytes)} bytes)"

        if magic_bytes != GZIP_MAGIC:
            return (
                f"File has .gz extension but does not have gzip magic bytes "
                f"(expected 0x1f 0x8b, found 0x{magic_bytes[0]:02x} 0x{magic_bytes[1]:02x})"
            )

    except OSError as e:
        return f"Failed to read file for gzip magic bytes check: {e}"

    return None


def run_grz_check(args: list[str]) -> subprocess.CompletedProcess:
    """
    Run `grz-check` with the given args.

    We catch KeyboardInterrupt here to allow the `grz-check` process to handle its own graceful shutdown
    and updating the progress logs from the report entries `grz-check` has generated so far.

    :param args: Arguments to pass to `grz-check`.
    :raises UserInterruptException: If interrupted via ctrl+c, `grz-check` exits with code 130.
    :raises subprocess.CalledProcessError: If `grz-check` fails for other reasons.
    """

    command = ["grz-check", *args]
    logging.info(f"Executing command: {' '.join(command)}")

    try:
        proc = subprocess.run(command, check=False)  # noqa: S603
    except KeyboardInterrupt:
        logging.warning("\nInterrupt received, allowing `grz-check` to shut down gracefully...")
        raise UserInterruptException from None

    if proc.returncode == 130:
        logging.warning("`grz-check` shut down gracefully due to an interrupt.")
        raise UserInterruptException
    elif proc.returncode != 0:
        logging.error(f"`grz-check` failed with a non-zero exit code: {proc.returncode}")
        raise subprocess.CalledProcessError(proc.returncode, command)
    return proc
