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
@grzcli.metadata_dir
@grzcli.logs_dir
@grzcli.encrypted_files_dir
@grzcli.threads
@grzcli.update_db
@click.option(
    "--consented/--non-consented",
    "consented",
    required=True,
    help="Whether to archive as consented or non-consented.",
)
def archive(  # noqa: PLR0913
    configuration: GrzctlConfig,
    submission_dir,
    metadata_dir,
    logs_dir,
    encrypted_files_dir,
    threads,
    update_db,
    consented,
    **kwargs,
):
    """
    Archive a submission within a GRZ/GDC.
    """
    bundled_mode = submission_dir is not None
    granular_mode = any(map(lambda v: v is not None, [metadata_dir, logs_dir, encrypted_files_dir]))

    if bundled_mode and granular_mode:
        raise click.UsageError("'--submission-dir' is mutually exclusive with explicit path options.")

    if bundled_mode:
        base = Path(submission_dir)
        _metadata_dir = base / "metadata"
        _logs_dir = base / "logs"
        _encrypted_files_dir = base / "encrypted_files"
    elif granular_mode:
        required = {
            "--metadata-dir": metadata_dir,
            "--logs-dir": logs_dir,
            "--encrypted-files-dir": encrypted_files_dir,
        }
        missing = [name for name, path in required.items() if path is None]
        if missing:
            raise click.UsageError(f"Granular mode requires: {', '.join(missing)}")
        _metadata_dir = Path(metadata_dir)
        _logs_dir = Path(logs_dir)
        _encrypted_files_dir = Path(encrypted_files_dir)
    else:
        raise click.UsageError("You must specify either '--submission-dir' or the required explicit path options.")

    archive_s3 = configuration.archives.consented.s3 if consented else configuration.archives.non_consented.s3

    log.info("Starting archival...")

    worker_inst = Worker(
        metadata_dir=_metadata_dir,
        files_dir=_metadata_dir.parent / "files",
        log_dir=_logs_dir,
        encrypted_files_dir=_encrypted_files_dir,
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
