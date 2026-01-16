"""Command for uploading a submission."""

import logging
from pathlib import Path

import click
from grz_common.utils.config import read_and_merge_config_files
from grz_common.workers.worker import Worker

log = logging.getLogger(__name__)

from grz_common.cli import config_file, read_config_from_ctx, submission_dir, threads

from ..models.config import UploadConfig


@click.command()
@submission_dir
@config_file
@threads
@click.pass_context
def upload(
    ctx,
    submission_dir,
    config_file: list[Path],
    threads,
):
    """
    Upload a submission to a GRZ/GDC.
    """
    config = UploadConfig.model_validate(read_config_from_ctx(ctx))

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
    submission_id = worker_inst.upload(config.s3)
    log.info(f"Generated submission ID for upload: {submission_id}")
    click.echo(submission_id)

    log.info("Upload finished!")
