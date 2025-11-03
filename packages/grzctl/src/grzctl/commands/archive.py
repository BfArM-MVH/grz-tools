"""Command for archiving a submission."""

import logging
from pathlib import Path

import click
from grz_common.cli import config_file, submission_dir, threads
from grz_common.workers.worker import Worker

from ..models.config import ArchiveConfig

log = logging.getLogger(__name__)


@click.command()
@submission_dir
@config_file
@threads
def archive(submission_dir, config_file, threads):
    """
    Archive a pre-staged submission directory from within a GRZ/GDC.
    This command expects a directory containing redacted metadata, redacted logs,
    and encrypted files.
    """
    config = ArchiveConfig.from_path(config_file)
    log.info("Starting archival upload...")

    submission_path = Path(submission_dir)

    worker_inst = Worker(
        metadata_dir=submission_path / "metadata",
        files_dir=submission_path / "files",
        log_dir=submission_path / "logs",
        encrypted_files_dir=submission_path / "encrypted_files",
        threads=threads,
    )
    worker_inst.archive(config.s3)

    log.info("Archival finished!")
