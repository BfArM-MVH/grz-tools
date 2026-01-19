"""Command for archiving a submission."""

import logging
from pathlib import Path
from typing import Any

import click
from grz_common.workers.worker import Worker

from ..models.config import ArchiveConfig

log = logging.getLogger(__name__)

import grz_common.cli as grzcli


@click.command()
@grzcli.configuration
@grzcli.submission_dir
@grzcli.threads
def archive(
    configuration: dict[str, Any],
    config_file: tuple[Path],
    submission_dir,
    threads,
):
    """
    Archive a submission within a GRZ/GDC.
    """
    config = ArchiveConfig.model_validate(configuration)

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
