"""Command for archiving a submission."""

import logging
from pathlib import Path
from typing import Any

import click
from grz_common.workers.worker import Worker
from grz_db.models.submission import SubmissionStateEnum

from ..dbcontext import DbContext
from ..models.config import ArchiveConfig

log = logging.getLogger(__name__)

import grz_common.cli as grzcli


@click.command()
@grzcli.configuration
@grzcli.submission_dir
@grzcli.threads
@grzcli.update_db
def archive(
    configuration: dict[str, Any],
    submission_dir,
    threads,
    update_db,
    **kwargs,
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
    submission_id = worker_inst.parse_encrypted_submission().submission_id
    with DbContext(
        configuration=configuration,
        submission_id=submission_id,
        start_state=SubmissionStateEnum.ARCHIVING,
        end_state=SubmissionStateEnum.ARCHIVED,
        enabled=update_db,
    ):
        worker_inst.archive(config.s3)

    log.info("Archival finished!")
