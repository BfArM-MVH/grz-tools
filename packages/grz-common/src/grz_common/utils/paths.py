"""Path utilities."""

import logging
import os
import stat
from os import PathLike
from pathlib import Path

log = logging.getLogger(__name__)


def is_relative_subdirectory(relative_path: str | PathLike, root_directory: str | PathLike) -> bool:
    """
    Check if the target path is a subdirectory of the root path
    using os.path.commonpath() without checking the file system.

    :param relative_path: The target path.
    :param root_directory: The root directory.
    :return: True if relative_path is a subdirectory of root_directory, otherwise False.
    """
    # Convert both paths to absolute paths without resolving symlinks
    root_directory = os.path.abspath(root_directory)
    relative_path = os.path.abspath(relative_path)

    common_path = os.path.commonpath([root_directory, relative_path])

    # Check if the common path is equal to the root path
    return common_path == root_directory


def ensure_directory_mode(base: Path, directory: Path, mode: int) -> None:
    """Create a directory below *base*, and set the permission bits of each directory on the way to *mode*.

    The walk goes from *base* down to *directory*.
    It creates each missing directory, and it runs ``chmod`` on each directory whose mode differs.
    ``chmod`` ignores the umask, so every directory ends up with exactly *mode*.
    *base* itself and the directories above it keep their mode.

    :param base: Existing directory that *directory* lies below. The helper does not touch it.
    :param directory: Directory to create and to set the mode of, together with its parents below *base*.
    :param mode: Permission bits, such as ``0o770``.
    :raises ValueError: If *directory* does not lie below *base*.
    """
    current = base
    for part in directory.relative_to(base).parts:
        current /= part
        current.mkdir(mode=mode, exist_ok=True)
        if stat.S_IMODE(current.stat().st_mode) != mode:
            current.chmod(mode)
