"""
Module: logging

This module provides functions for setting up logging configuration.
"""

from __future__ import annotations

import logging
import sys
from os import PathLike
from pathlib import Path

try:
    from tqdm.auto import tqdm

    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False

log = logging.getLogger(__name__)

LOGGING_FORMAT = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
LOGGING_DATEFMT = "%Y-%m-%d %I:%M %p"


class TqdmLoggingHandler(logging.Handler):
    """
    A logging handler that outputs to stderr via tqdm.write().
    This ensures log messages don't interfere with tqdm progress bars.
    """

    def emit(self, record):
        try:
            msg = self.format(record)
            tqdm.write(msg, file=sys.stderr)
            self.flush()
        except Exception:
            self.handleError(record)


def add_filelogger(file_path: str | PathLike, level: str = "INFO", logger_name: str | None = None) -> None:
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
                        Default is the root logger.
    """
    # passing None to getLogger gets the root logger
    logger = logging.getLogger(logger_name)

    file_path = Path(file_path)

    file_handler = logging.FileHandler(file_path)
    file_handler.setLevel(level.upper())
    file_handler.setFormatter(logging.Formatter(fmt=LOGGING_FORMAT, datefmt=LOGGING_DATEFMT))
    logger.addHandler(file_handler)
    log.info(
        "File logger added for %s at %s with level %s.",
        logger.name,
        file_path,
        level.upper(),
    )


def setup_cli_logging(log_file: str | None, log_level: str):
    """
    Setup logging for the CLI.
    Uses TqdmLoggingHandler if available to play nicely with progress bars.
    """
    root_logger = logging.getLogger()
    root_logger.setLevel(log_level.upper())

    formatter = logging.Formatter(fmt=LOGGING_FORMAT, datefmt=LOGGING_DATEFMT)

    # Remove existing handlers to avoid duplication (e.g. default handlers)
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)

    console_handler = TqdmLoggingHandler() if HAS_TQDM else logging.StreamHandler(sys.stderr)

    console_handler.setFormatter(formatter)
    root_logger.addHandler(console_handler)

    if log_file:
        # add file handler to root logger
        add_filelogger(
            log_file,
            log_level.upper(),
        )

    log.debug("Logging setup complete.")
