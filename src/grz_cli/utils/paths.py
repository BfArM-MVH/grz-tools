"""Path utilities."""

import logging
import os
from os import PathLike

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
