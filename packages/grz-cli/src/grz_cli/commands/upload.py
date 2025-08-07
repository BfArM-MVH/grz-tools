"""Command for uploading a submission."""

import logging
from pathlib import Path

import click
from grz_common.workers.upload.boto import S3BotoUploadWorker
from grz_common.workers.upload.gatekeeper import GrzGatekeeperUploadWorker
from grz_common.workers.worker import Worker

log = logging.getLogger(__name__)

from grz_common.cli import config_file, submission_dir, threads

from ..models.config import UploadConfig


@click.command()
@submission_dir
@config_file
@threads
def upload(submission_dir, config_file, threads):
    """
    Upload a submission to a GRZ/GDC.

    Uses the grz-gatekeeper API if an 'api_base_url' is configured,
    otherwise falls back to the direct S3 Boto worker.
    """
    config = UploadConfig.from_path(config_file)

    log.info("Starting upload...")

    submission_dir = Path(submission_dir)
    log_dir = submission_dir / "logs"

    if api_base_url := config.s3.api_base_url:
        log.info(f"API endpoint configured. Using GrzGatekeeperUploadWorker with {api_base_url}")
        upload_worker = GrzGatekeeperUploadWorker(
            api_base_url=api_base_url,
            status_file_path=log_dir / "progress_upload.cjson",  # or use different file?
            threads=threads,
        )
    else:
        log.info("No API endpoint configured. Falling back to direct S3BotoUploadWorker.")
        upload_worker = S3BotoUploadWorker(
            s3_options=config.s3,
            status_file_path=log_dir / "progress_upload.cjson",
            threads=threads,
        )

    worker = Worker(
        metadata_dir=submission_dir / "metadata",
        files_dir=submission_dir / "files",
        log_dir=log_dir,
        encrypted_files_dir=submission_dir / "encrypted_files",
    )
    # output the generated submission ID
    submission_id = worker.upload(upload_worker=upload_worker)
    log.info(f"Generated submission ID for upload: {submission_id}")
    click.echo(submission_id)

    log.info("Upload finished!")
