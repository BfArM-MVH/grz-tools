"""Command for uploading a submission."""

import logging
from pathlib import Path

import click
from grz_common.workers.worker import Worker

log = logging.getLogger(__name__)

from grz_common.cli import config_file, submission_dir, threads

from ..models.config import UploadConfig


@click.command()
@submission_dir
@config_file
@threads
def upload(
    submission_dir,
    config_file,
    threads,
):
    """
    Upload a submission to a GRZ/GDC.
    """
    config = UploadConfig.from_path(config_file)

    log.info("Starting upload...")

    submission_dir = Path(submission_dir)

    worker_inst = Worker(
        metadata_dir=submission_dir / "metadata",
        files_dir=submission_dir / "files",
        log_dir=submission_dir / "logs",
        encrypted_files_dir=submission_dir / "encrypted_files",
        threads=threads,
    )
    # output the generated submission ID
    click.echo(worker_inst.upload(config.s3))

    log.info("Upload finished!")
