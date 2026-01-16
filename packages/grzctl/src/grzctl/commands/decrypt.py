"""Command for decrypting a submission."""

import logging
import sys
from pathlib import Path

import click
from grz_common.cli import config_file, config_files_from_ctx, force, submission_dir
from grz_common.utils.config import read_and_merge_config_files
from grz_common.workers.worker import Worker

from ..models.config import DecryptConfig

log = logging.getLogger(__name__)


@click.command()
@submission_dir
@config_file
@force
@click.pass_context
def decrypt(ctx, submission_dir, config_file: list[Path], force):
    """
    Decrypt a submission.

    Decrypting a submission requires the _private_ key of the original recipient.
    """
    config_files = config_files_from_ctx(ctx)
    config = DecryptConfig.model_validate(read_and_merge_config_files(config_files))

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
