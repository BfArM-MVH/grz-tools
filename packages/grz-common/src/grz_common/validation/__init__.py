"""Validation package for grz-cli."""

import subprocess
from importlib.resources import as_file, files


def run_grz_check(args: list[str]) -> subprocess.CompletedProcess:
    """
    Finds and executes the 'grz-check' binary with the given arguments.

    Returns:
        subprocess.CompletedProcess: The result of the execution.

    Raises:
        FileNotFoundError: If the 'grz-check' binary cannot be found.
        subprocess.CalledProcessError: If the binary returns a non-zero exit code.
    """
    binary_name = "grz-check"

    try:
        binary_path = files("grz_check_provider.bin").joinpath(binary_name)
    except ModuleNotFoundError as e:
        raise FileNotFoundError("Could not find the 'grz-check-provider' package. Is it installed?") from e

    with as_file(binary_path) as binary_path:
        command = [str(binary_path), *args]
        print(f"Executing command: {' '.join(command)}")

        binary_path.chmod(0o755)

        return subprocess.run(
            " ".join(command),
            capture_output=True,
            text=True,
            check=True,
            shell=True,
        )
