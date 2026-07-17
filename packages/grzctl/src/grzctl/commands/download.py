"""Command for downloading a submission."""

import logging
from pathlib import Path

import click
import grz_common.cli as grzcli
from grz_common.transfer import get_metadata_upload_timestamp, init_s3_client
from grz_common.workers.worker import Worker
from grz_db.models.submission import SubmissionStateEnum

from ..commands import grzctl_configuration
from ..dbcontext import DbContext
from ..models.config import GrzctlConfig

log = logging.getLogger(__name__)


@click.command()
@grzctl_configuration
@grzcli.submission_id
@grzcli.output_dir
@grzcli.threads
@grzcli.force
@grzcli.update_db
@click.option(
    "--inbox-bucket",
    default=None,
    help="Inbox bucket name to use. Required when a submitter has multiple inboxes configured.",
)
@click.option(
    "--populate/--no-populate",
    default=True,
    help="Update the submission metadata with information from metadata.json and S3. If combined with --force, will overwrite information in db without asking.",
)
def download(  # noqa: PLR0913
    configuration: GrzctlConfig,
    submission_id,
    output_dir,
    threads,
    force,
    update_db,
    inbox_bucket,
    populate,
    **kwargs,
):
    """
    Download a submission from a GRZ.

    Downloaded metadata is stored within the `metadata` sub-folder of the submission output directory.
    Downloaded files are stored within the `encrypted_files` sub-folder of the submission output directory.
    """
    config = configuration
    s3_options = config.resolve_inbox_by_submission_id(submission_id, inbox_bucket).s3

    log.info("Starting download...")

    submission_dir_path = Path(output_dir)
    if not submission_dir_path.is_dir():
        log.debug("Creating submission directory %s", submission_dir_path)
        submission_dir_path.mkdir(mode=0o770, parents=False, exist_ok=False)

    worker_inst = Worker(
        metadata_dir=submission_dir_path / "metadata",
        files_dir=submission_dir_path / "files",
        log_dir=submission_dir_path / "logs",
        encrypted_files_dir=submission_dir_path / "encrypted_files",
        threads=threads,
    )

    with DbContext(
        configuration=configuration,
        submission_id=submission_id,
        start_state=SubmissionStateEnum.DOWNLOADING,
        end_state=SubmissionStateEnum.DOWNLOADED,
        enabled=update_db,
    ) as db_context:
        worker_inst.download(s3_options, submission_id, force=force)
        if populate:
            if not db_context.db:
                log.warning("Database context is not available, skipping population of submission metadata in DB.")
            else:
                s3_client = init_s3_client(s3_options)
                submission_date = get_metadata_upload_timestamp(s3_client, s3_options.bucket, submission_id).date()
                metadata = worker_inst.parse_submission().metadata.content
                db_context.db.populate(
                    submission_id,
                    metadata,
                    submission_date,
                    force=force,
                    on_missing="create",
                )

    log.info("Download finished!")
