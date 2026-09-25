"""Command for downloading a submission."""

import logging
from pathlib import Path

import click
import grz_common.cli as grzcli
from grz_common.transfer import get_metadata_upload_timestamp, init_s3_client
from grz_common.utils.version_check import check_metadata_version_and_exit_if_needed
from grz_common.workers.worker import Worker
from grz_db.errors import SubmissionNotFoundError
from grz_db.models.submission import SubmissionStateEnum

from ..commands import grzctl_configuration, inbox_option
from ..dbcontext import DbContext
from ..models.config import GrzctlConfig
from .db.cli import get_submission_db_or_none
from .inbox_resolution import require_inbox

log = logging.getLogger(__name__)


@click.command()
@grzctl_configuration
@grzcli.submission_id
@grzcli.output_dir
@grzcli.threads
@grzcli.force
@grzcli.update_db
@inbox_option
@click.option(
    "--populate/--no-populate",
    default=True,
    help="Update the submission metadata with information from metadata.json and S3.",
)
def download(  # noqa: PLR0913, PLR0917
    configuration: GrzctlConfig,
    submission_id: str,
    output_dir: str,
    threads: int,
    force: bool,
    update_db: bool,
    inbox_name: str,
    populate: bool,
    **kwargs,
):
    """
    Download a submission from a GRZ.

    Downloaded metadata is stored within the `metadata` sub-folder of the submission output directory.
    Downloaded files are stored within the `encrypted_files` sub-folder of the submission output directory.
    """
    submitter_id = submission_id.split("_", maxsplit=1)[0]
    resolved_inbox = require_inbox(
        configuration,
        submitter_id=submitter_id,
        submission_id=submission_id,
        inbox_name=inbox_name,
        db_service=get_submission_db_or_none(configuration),
        scan=True,
    )
    s3_options = configuration.resolve_inbox(submitter_id=submitter_id, inbox_name=resolved_inbox).s3
    bucket_name = s3_options.bucket
    inbox_desc = f"'{resolved_inbox}' (bucket '{bucket_name}')" if resolved_inbox != bucket_name else f"'{bucket_name}'"

    log.info(f"Starting download from inbox {inbox_desc}...")

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
        worker_inst.download(
            s3_options,
            submission_id,
            force=force,
            metadata_version_check=lambda metadata_schema_version: check_metadata_version_and_exit_if_needed(
                s3_options,
                metadata_schema_version,
            ),
        )
        if db_context.db:
            if populate:
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
            # The download knows the inbox, so record it.
            # Commands and the decryption key lookup can then find the submission without
            # being told the inbox again.
            log.info(f"Recording inbox {resolved_inbox} for submission {submission_id}...")
            try:
                db_context.db.set_submission_inbox(submission_id, resolved_inbox)
            except SubmissionNotFoundError:
                # Without populate, a download does not register the submission in the database.
                # There is then no row to record the inbox on.
                # A later populate or sync-from-inbox will record it.
                log.debug("Submission %s is not registered in the database, skipping inbox recording.", submission_id)
        elif populate:
            log.warning("Database context is not available, skipping population of submission metadata in DB.")

    log.info("Download finished!")
