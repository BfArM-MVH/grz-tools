"""
Module: logging_setup

This module provides functions for setting up logging configuration.
"""

from __future__ import annotations

import logging
from os import PathLike
from pathlib import Path

from .constants import LOGGING_DATEFMT, LOGGING_FORMAT, PACKAGE_ROOT

log = logging.getLogger(__name__)


def add_filelogger(
    file_path: str | PathLike | None = None,
    level: str = "INFO",
    logger_name: str = PACKAGE_ROOT,
) -> None:
    """
    Add file logging for the specified package.

    This function configures a file logger to capture log messages
    for the package specified by logger_name. If no file path
    is provided, a default log file will be created in the user's
    home directory.

    :param file_path: Optional; the path to the log file. If None,
                      a default path will be used.
    :param level: Optional; the logging level. Default is 'INFO'.
                  Must be a valid logging level name (e.g., 'DEBUG', 'INFO').
    :param logger_name: Optional; the name of the logger to add the file handler to.
                        Default is PACKAGE_ROOT.
    """
    logger = logging.getLogger(logger_name)

    if file_path is None:
        default_log_dir = Path.home() / "logs"
        default_log_dir.mkdir(parents=True, exist_ok=True)
        file_path = default_log_dir / f"{logger_name}.log"
        log.warning("No log file path provided, using default: %s", file_path)
    else:
        file_path = Path(file_path)

    try:
        fh = logging.FileHandler(file_path)
        fh.setLevel(level.upper())
        fh.setFormatter(logging.Formatter(LOGGING_FORMAT, LOGGING_DATEFMT))
        logger.addHandler(fh)
        log.info(
            "File logger added for %s at %s with level %s.",
            logger.name,
            file_path,
            level.upper(),
        )
    except Exception as e:
        log.error("Failed to add file logger: %s", e)
