"""Command for archiving a submission."""

import logging
from pathlib import Path

import click
import grz_common.cli as grzcli
from grz_common.workers.worker import Worker
from grz_db.models.submission import SubmissionStateEnum

from ..commands import grzctl_configuration
from ..dbcontext import DbContext
from ..models.config import GrzctlConfig

log = logging.getLogger(__name__)


@click.command()
@grzctl_configuration
@grzcli.submission_dir
@grzcli.threads
@grzcli.update_db
@click.option(
    "--consented/--non-consented",
    "consented",
    default=True,
    help="Whether to archive as consented (default) or non-consented.",
)
def archive(
    configuration: GrzctlConfig,
    submission_dir,
    threads,
    update_db,
    consented,
    **kwargs,
):
    """
    Archive a submission within a GRZ/GDC.
    """
    config = configuration

    archive_s3 = config.archives.consented.s3 if consented else config.archives.non_consented.s3

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
        worker_inst.archive(archive_s3)

    log.info("Archival finished!")
