"""Command for decrypting a submission."""

import logging
import sys
from pathlib import Path
from typing import Any

import click
import grz_common.cli as grzcli
from grz_common.workers.worker import Worker
from grz_db.models.submission import SubmissionStateEnum

from ..dbcontext import DbContext
from ..models.config import DecryptConfig

log = logging.getLogger(__name__)


@click.command()
@grzcli.configuration
@grzcli.submission_dir
@grzcli.force
@grzcli.update_db
def decrypt(
    configuration: dict[str, Any],
    config_file: tuple[Path],
    submission_dir,
    force,
    update_db,
    **kwargs,
):
    """
    Decrypt a submission.

    Decrypting a submission requires the _private_ key of the original recipient.
    """
    config = DecryptConfig.model_validate(configuration)

    grz_privkey_path = config.keys.grz_private_key_path
    if not grz_privkey_path:
        log.error("GRZ private key path is required for decryption.")
        sys.exit(1)

    log.info("Starting decryption...")

    submission_dir = Path(submission_dir)

    worker_inst = Worker(
        metadata_dir=submission_dir / "metadata",
        files_dir=submission_dir / "files",
        log_dir=submission_dir / "logs",
        encrypted_files_dir=submission_dir / "encrypted_files",
    )
    submission_id = worker_inst.parse_encrypted_submission().submission_id
    with DbContext(
        configuration=configuration,
        submission_id=submission_id,
        start_state=SubmissionStateEnum.DECRYPTING,
        end_state=SubmissionStateEnum.DECRYPTED,
        enabled=update_db,
    ):
        worker_inst.decrypt(grz_privkey_path, force=force)

    log.info("Decryption successful!")
