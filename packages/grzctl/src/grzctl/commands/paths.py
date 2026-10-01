"""Resolution of submission directories from a bundled or explicit options.

Commands that touch one part of a submission accept either one bundled directory
(``--submission-dir`` / ``--output-dir``) or the individual directory options, so that a
caller can keep parts of a submission in directories of its own choosing. This module
resolves whichever form the caller used into one mapping, so that a command states which
directories it needs instead of repeating the mode handling.
"""

from pathlib import Path
from typing import Any

import click

SUBDIR_BY_OPTION: dict[str, str] = {
    "--metadata-dir": "metadata",
    "--files-dir": "files",
    "--encrypted-files-dir": "encrypted_files",
    "--logs-dir": "logs",
    "--output-encrypted-files-dir": "encrypted_files",
}
"""The subdirectory of a bundled directory that each explicit option names."""


def resolve_dirs(bundled_dir: Any, bundled_option: str, explicit: dict[str, Any]) -> dict[str, Path]:
    """
    Resolve a bundled directory or a set of explicit directory options into one mapping.

    Every name in *explicit* is resolved, so the returned mapping has the same keys. A
    command asks for exactly the directories it needs, and this function decides whether
    they come from one bundled directory or from the individual options.

    :param bundled_dir: Value of the bundled directory option, or ``None`` if it was not given.
    :param bundled_option: Name of that option, for error messages.
    :param explicit: Explicit directory options, keyed by option name as declared in
        :data:`SUBDIR_BY_OPTION`.
    :returns: Mapping of option name to resolved path.
    :raises click.UsageError: If both forms are given, if neither is, or if an explicit
        option is missing.
    """
    given = {name: value for name, value in explicit.items() if value is not None}

    if bundled_dir is not None:
        if given:
            raise click.UsageError(f"'{bundled_option}' is mutually exclusive with the explicit path options.")
        base = Path(bundled_dir)
        return {name: base / SUBDIR_BY_OPTION[name] for name in explicit}

    missing = [name for name in explicit if name not in given]
    if missing:
        raise click.UsageError(
            f"You must specify either '{bundled_option}' or the explicit path options: {', '.join(missing)}."
        )
    return {name: Path(value) for name, value in explicit.items()}
