"""Module providing functions for validating submissions, including parallel checksum calculation across files"""

import logging
import sys
from collections.abc import Generator
from os import PathLike
from pathlib import Path

from tqdm.contrib.concurrent import process_map

from ..progress_logging import FileProgressLogger
from ..submission import SubmissionFileMetadata

type ProcessItem = tuple[Path, SubmissionFileMetadata]
type ProcessResult = tuple[Path, SubmissionFileMetadata, dict]


def validate_checksums(
    files, progress_log_file: str | PathLike, logger: logging.Logger
) -> Generator[str, None, None]:
    """
    Validates the checksum of the files against the metadata and prints the errors.

    :return: Generator of errors
    """
    progress_logger = FileProgressLogger(log_file_path=progress_log_file)
    # cleanup log file and keep only files listed here
    progress_logger.cleanup(
        keep=[(file_path, file_metadata) for file_path, file_metadata in files.items()]
    )
    yield from _parallel_validate(files, progress_logger, logger)


def _parallel_validate(files, progress_logger, logger) -> list[str]:
    files_to_validate = _determine_files_to_validate(files, progress_logger, logger)

    for file_path, file_metadata, state in process_map(
        _validate_item, enumerate(files_to_validate)
    ):
        progress_logger.set_state(file_path, file_metadata, state=state)
        if state["validation_passed"] is False:
            # Do we want to exit early or collect all errors and then exit?
            sys.exit(f"Validation failed for {file_path}, exiting.")
    # If we exit early, there's nothing to report here, really.
    return []


def _determine_files_to_validate(
    files, progress_logger, logger: logging.Logger
) -> list[ProcessItem]:
    logged_errors = []
    files_to_validate = []
    for local_file_path, file_metadata in files.items():
        logged_state = progress_logger.get_state(local_file_path, file_metadata)

        # determine if we can skip the verification
        if logged_state is None:
            logger.debug("State for %s not calculated yet", local_file_path)
        elif not logged_state.get("validation_passed", False):
            errors = logged_state.get("errors", [])
            logged_errors.extend(errors)

            # skip re-verification
            continue
        else:
            logger.debug(
                "Validation for %s already passed, skipping...",
                str(local_file_path),
            )

            # skip re-verification
            continue

        logger.debug("Validating '%s'...", str(local_file_path))
        # validate the file
        files_to_validate.append((local_file_path, file_metadata))
    if logged_errors:
        logger.error("Validation errors found:")
        for error in logged_errors:
            logger.error(error)
        sys.exit("Validation errors found, exiting.")

    return files_to_validate


def _validate_item(item: tuple[int, ProcessItem]) -> ProcessResult:
    idx, (file_path, file_metadata) = item
    errors = list(
        file_metadata.validate_data(file_path, tqdm_kwargs=dict(position=idx))
    )
    validation_passed = len(errors) == 0
    return (
        file_path,
        file_metadata,
        {
            "errors": errors,
            "validation_passed": validation_passed,
        },
    )
