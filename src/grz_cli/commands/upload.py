"""Command for uploading a submission."""

import logging
from pathlib import Path

import click

from ..workers.worker import Worker

log = logging.getLogger(__name__)

from ..utils.config import read_config
from .common import config_file, submission_dir, threads


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
    config = read_config(config_file)

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
    print(worker_inst.upload(config))

    log.info("Upload finished!")
