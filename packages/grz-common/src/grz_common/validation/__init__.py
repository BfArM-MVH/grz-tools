"""Validation package for grz-cli."""

import subprocess


def run_grz_check(args: list[str]) -> subprocess.CompletedProcess:
    """
    Finds and executes the 'grz-check' binary with the given arguments.

    Returns:
        subprocess.CompletedProcess: The result of the execution.

    Raises:
        FileNotFoundError: If the 'grz-check' binary cannot be found.
        subprocess.CalledProcessError: If the binary returns a non-zero exit code.
    """

    command = ["grz-check", *args]
    print(f"Executing command: {' '.join(command)}")

    return subprocess.run(  # noqa: S603
        command,
        capture_output=True,
        text=True,
        check=True,
    )
