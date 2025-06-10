"""Command for archiving a submission."""

import logging
from pathlib import Path

import click
from grz_common.workers.worker import Worker

from ..models.config import ArchiveConfig

log = logging.getLogger(__name__)

from grz_common.cli import config_file, submission_dir, threads


@click.command()
@submission_dir
@config_file
@threads
def archive(
    submission_dir,
    config_file,
    threads,
):
    """
    Archive a submission within a GRZ/GDC.
    """
    config = ArchiveConfig.from_path(config_file)

    log.info("Starting archival...")

    submission_dir = Path(submission_dir)

    worker_inst = Worker(
        metadata_dir=submission_dir / "metadata",
        files_dir=submission_dir / "files",
        log_dir=submission_dir / "logs",
        encrypted_files_dir=submission_dir / "encrypted_files",
        threads=threads,
    )
    worker_inst.archive(config.s3)

    log.info("Archival finished!")
