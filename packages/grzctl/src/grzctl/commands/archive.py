"""Command for archiving a submission."""

import logging

import click
import grz_common.cli as grzcli
from grz_common.workers.worker import Worker
from grz_db.models.submission import SubmissionStateEnum

from ..commands import grzctl_configuration
from ..dbcontext import DbContext
from ..models.config import GrzctlConfig
from .paths import encrypted_files_dir_option, logs_dir_option, metadata_dir_option, resolve_dirs, submission_dir_option

log = logging.getLogger(__name__)


@click.command()
@grzctl_configuration
@submission_dir_option
@metadata_dir_option
@logs_dir_option
@encrypted_files_dir_option
@grzcli.threads
@grzcli.update_db
def archive(  # noqa: PLR0913, PLR0917
    configuration: GrzctlConfig,
    submission_dir,
    metadata_dir,
    logs_dir,
    encrypted_files_dir,
    threads,
    update_db,
    **kwargs,
):
    """
    Archive a submission within a GRZ/GDC.
    """
    paths = resolve_dirs(
        bundled_dir=submission_dir,
        bundled_option="--submission-dir",
        explicit={
            "--metadata-dir": metadata_dir,
            "--logs-dir": logs_dir,
            "--encrypted-files-dir": encrypted_files_dir,
        },
    )
    metadata_path = paths["--metadata-dir"]

    log.info("Starting archival...")

    worker_inst = Worker(
        metadata_dir=metadata_path,
        files_dir=metadata_path.parent / "files",
        log_dir=paths["--logs-dir"],
        encrypted_files_dir=paths["--encrypted-files-dir"],
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
