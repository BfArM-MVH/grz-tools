"""Validation package for grz-cli."""
import logging
import subprocess


def run_grz_check(args: list[str]) -> subprocess.CompletedProcess:
    """
    Run `grz-check` with the given args.

    :param args: Arguments to pass to `grz-check`.
    """

    command = ["grz-check", *args]
    logging.info(f"Executing command: {' '.join(command)}")

    return subprocess.run(  # noqa: S603
        command,
        capture_output=False,
        text=True,
        check=True,
    )
