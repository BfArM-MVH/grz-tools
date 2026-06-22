"""
Module: logging

This module provides functions for setting up logging configuration.
"""

from __future__ import annotations

import logging
import socket
from os import PathLike
from pathlib import Path

log = logging.getLogger(__name__)


def get_hostname() -> str | None:
    """
    Return the current host name, or ``None`` if it cannot be determined.

    :return: The host name, or ``None`` if :func:`socket.gethostname` fails.
    """
    try:
        return socket.gethostname()
    except OSError:
        return None


def build_logging_format(hostname: str | None) -> str:
    """
    Build the logging format string, including the host name when available.

    :param hostname: The host name to embed, or ``None`` to omit it.
    :return: A :mod:`logging` format string.
    """
    hostname_part = f"{hostname} " if hostname is not None else ""
    return f"%(asctime)s [%(levelname)s] {hostname_part}%(name)s: %(message)s"


HOSTNAME: str | None = get_hostname()
LOGGING_FORMAT = build_logging_format(HOSTNAME)
LOGGING_DATEFMT = "%Y-%m-%d %I:%M %p"


class AlembicInfoNoiseFilter(logging.Filter):
    """Filter out repetitive initialization messages from Alembic."""

    def filter(self, record: logging.LogRecord) -> bool:
        if record.name == "alembic.runtime.migration":
            msg = record.getMessage()
            if "Context impl" in msg or "Will assume transactional DDL" in msg:
                return False
        return True


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
    # set the root log level since this is the CLI
    logging.getLogger().setLevel(log_level.upper())

    logging.basicConfig(level=log_level.upper(), format=LOGGING_FORMAT, datefmt=LOGGING_DATEFMT)

    if log_file:
        # add file handler to root logger
        add_filelogger(
            log_file,
            log_level.upper(),
        )
    logging.getLogger("alembic.runtime.migration").addFilter(AlembicInfoNoiseFilter())

    log.debug("Logging setup complete.")
