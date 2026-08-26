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
def archive(
    configuration: GrzctlConfig,
    submission_dir,
    threads,
    update_db,
    **kwargs,
):
    """
    Archive a submission within a GRZ/GDC.
    """
    log.info("Starting archival...")

    submission_dir = Path(submission_dir)

    worker_inst = Worker(
        metadata_dir=submission_dir / "metadata",
        files_dir=submission_dir / "files",
        log_dir=submission_dir / "logs",
        encrypted_files_dir=submission_dir / "encrypted_files",
        threads=threads,
    )
    encrypted_submission = worker_inst.parse_encrypted_submission()
    submission_id = encrypted_submission.submission_id

    submission_date = encrypted_submission.metadata.content.submission.submission_date
    consented = encrypted_submission.metadata.content.consents_to_research(submission_date)

    archive_s3 = configuration.archives.consented.s3 if consented else configuration.archives.non_consented.s3

    with DbContext(
        configuration=configuration,
        submission_id=submission_id,
        start_state=SubmissionStateEnum.ARCHIVING,
        end_state=SubmissionStateEnum.ARCHIVED,
        enabled=update_db,
    ):
        worker_inst.archive(archive_s3)

    log.info("Archival finished!")
