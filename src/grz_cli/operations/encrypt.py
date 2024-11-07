"""Module providing functions for encrypting submissions, including parallel encryption across files"""

import logging
import sys
from collections.abc import Generator
from functools import partial
from pathlib import Path

from tqdm.contrib.concurrent import process_map

from ..file_operations import Crypt4GH
from ..progress_logging import FileProgressLogger
from ..submission import EncryptedSubmission, Submission, SubmissionFileMetadata

type ProcessItem = tuple[Path, Path, SubmissionFileMetadata, dict, logging.Logger]
type ProcessResult = tuple[Path, SubmissionFileMetadata, dict]


def _parallel_encrypt(  # noqa: PLR0913
    submission: Submission,
    public_keys,
    encrypted_files_dir: Path,
    progress_logger: FileProgressLogger,
    threads: int | None = None,
    logger: logging.Logger = logging.getLogger(__name__),
):
    files_to_encrypt = list(
        _determine_files_to_encrypt(
            submission, encrypted_files_dir, progress_logger, logger
        )
    )

    for file_path, file_metadata, state in process_map(
        partial(_encrypt_item, public_keys),
        enumerate(files_to_encrypt),
        max_workers=threads,
    ):
        progress_logger.set_state(file_path, file_metadata, state=state)
        if state["encryption_successful"] is False:
            sys.exit(f"Encryption failed for {file_path}, exiting.")


def _determine_files_to_encrypt(
    submission: Submission, encrypted_files_dir, progress_logger, logger: logging.Logger
) -> Generator[ProcessItem, None, None]:
    for file_path, file_metadata in submission.files.items():
        logged_state = progress_logger.get_state(file_path, file_metadata)
        logger.debug("state for %s: %s", file_path, logged_state)

        encrypted_file_path = (
            encrypted_files_dir
            / EncryptedSubmission.get_encrypted_file_path(file_metadata.file_path)
        )
        encrypted_file_path.parent.mkdir(mode=0o770, parents=True, exist_ok=True)

        if (
            (logged_state is None)
            or not logged_state.get("encryption_successful", False)
            or not encrypted_file_path.is_file()
        ):
            yield (
                file_path,
                encrypted_file_path,
                file_metadata,
                logged_state,
                logger,
            )
        else:
            logger.info(
                "File '%s' already encrypted in '%s'",
                str(file_path),
                str(encrypted_file_path),
            )


def _encrypt_item(public_keys, item: tuple[int, ProcessItem]) -> ProcessResult:
    idx, (file_path, encrypted_file_path, file_metadata, logged_state, logger) = item
    logger.debug("state for %s: %s", file_path, logged_state)

    logger.info(
        "Encrypting file: '%s' -> '%s'",
        str(file_path),
        str(encrypted_file_path),
    )

    try:
        Crypt4GH.encrypt_file(
            file_path,
            encrypted_file_path,
            public_keys,
            tqdm_kwargs=dict(position=idx, leave=True),
        )

        logger.info(f"Encryption complete for {str(file_path)}. ")
        return file_path, file_metadata, {"encryption_successful": True}
    except Exception as e:
        logger.error("Encryption failed for '%s'", str(file_path))

        return (
            file_path,
            file_metadata,
            {"encryption_successful": False, "error": str(e)},
        )
