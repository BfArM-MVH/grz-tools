"""Command for downloading a submission."""

import logging
from pathlib import Path

import click
import grz_common.cli as grzcli
from grz_cli.utils.version_check import check_metadata_version_and_exit_if_needed
from grz_common.transfer import get_metadata_upload_timestamp, init_s3_client
from grz_common.workers.worker import Worker
from grz_db.models.submission import SubmissionStateEnum

from ..commands import grzctl_configuration, inbox_option
from ..dbcontext import DbContext
from ..models.config import GrzctlConfig

log = logging.getLogger(__name__)


@click.command()
@grzctl_configuration
@grzcli.submission_id
@click.option(
    "--output-dir",
    "output_dir",
    type=grzcli.DIR_RW_C,
    default=None,
    help="Path to the target submission output directory.",
)
@grzcli.metadata_dir
@grzcli.encrypted_files_dir
@grzcli.logs_dir
@grzcli.threads
@grzcli.force
@grzcli.update_db
@inbox_option
@click.option(
    "--populate/--no-populate",
    default=True,
    help="Update the submission metadata with information from metadata.json and S3.",
)
def download(  # noqa: PLR0913
    configuration: GrzctlConfig,
    submission_id: str,
    output_dir: str,
    metadata_dir,
    encrypted_files_dir,
    logs_dir,
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
    bundled_mode = output_dir is not None
    granular_mode = any(map(lambda v: v is not None, [metadata_dir, encrypted_files_dir, logs_dir]))

    if bundled_mode and granular_mode:
        raise click.UsageError("'--output-dir' is mutually exclusive with explicit path options.")

    if bundled_mode:
        base = Path(output_dir)
        _metadata_dir = base / "metadata"
        _encrypted_files_dir = base / "encrypted_files"
        _logs_dir = base / "logs"
    elif granular_mode:
        required = {
            "--metadata-dir": metadata_dir,
            "--encrypted-files-dir": encrypted_files_dir,
            "--logs-dir": logs_dir,
        }
        missing = [name for name, path in required.items() if path is None]
        if missing:
            raise click.UsageError(f"Granular mode requires: {', '.join(missing)}")
        _metadata_dir = Path(metadata_dir)
        _encrypted_files_dir = Path(encrypted_files_dir)
        _logs_dir = Path(logs_dir)
    else:
        raise click.UsageError("You must specify either '--output-dir' or the required explicit path options.")

    submitter_id = submission_id.split("_", maxsplit=1)[0]
    s3_options = configuration.resolve_inbox(submitter_id=submitter_id, inbox_name=inbox_name).s3
    bucket_name = s3_options.bucket
    inbox_desc = f"'{inbox_name}' (bucket '{bucket_name}')" if inbox_name != bucket_name else f"'{bucket_name}'"

    log.info(f"Starting download from inbox {inbox_desc}...")

    if bundled_mode:
        submission_dir_path = base
        if not submission_dir_path.is_dir():
            log.debug("Creating submission directory %s", submission_dir_path)
            submission_dir_path.mkdir(mode=0o770, parents=False, exist_ok=False)

    worker_inst = Worker(
        metadata_dir=_metadata_dir,
        files_dir=_metadata_dir.parent / "files",
        log_dir=_logs_dir,
        encrypted_files_dir=_encrypted_files_dir,
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
