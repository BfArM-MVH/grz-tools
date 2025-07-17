"""Validation package for grz-cli."""

import logging
import subprocess


class UserInterruptException(Exception):
    """Raised when an interrupt is triggered via ctrl+c."""

    pass


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
