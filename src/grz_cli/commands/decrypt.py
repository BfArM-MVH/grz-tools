"""Command for decrypting a submission."""

import logging
import sys
from pathlib import Path

import click

from ..utils.config import read_config
from ..workers.worker import Worker
from .common import config_file, submission_dir

log = logging.getLogger(__name__)


@click.command()
@submission_dir
@config_file
def decrypt(
    submission_dir,
    config_file,
):
    """
    Decrypt a submission.

    Decrypting a submission requires the _private_ key of the original recipient.
    """
    config = read_config(config_file)

    grz_privkey_path = config.grz_private_key_path
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
    worker_inst.decrypt(grz_privkey_path)

    log.info("Decryption successful!")
