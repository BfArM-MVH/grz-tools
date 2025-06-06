"""Command for decrypting a submission."""

import logging
import sys
from pathlib import Path

import click
from grz_common.cli import config_file, force, submission_dir
from grz_common.workers.worker import Worker

from ..models.config import DecryptConfig

log = logging.getLogger(__name__)


@click.command()
@submission_dir
@config_file
@force
def decrypt(submission_dir, config_file, force):
    """
    Decrypt a submission.

    Decrypting a submission requires the _private_ key of the original recipient.
    """
    config = DecryptConfig.from_path(config_file)

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
    worker_inst.decrypt(grz_privkey_path, force=force)

    log.info("Decryption successful!")
