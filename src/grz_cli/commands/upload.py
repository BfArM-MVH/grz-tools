"""Command for uploading a submission."""

import logging
from os import sched_getaffinity
from pathlib import Path

import click

from ..workers.worker import Worker

log = logging.getLogger(__name__)


@click.command()
@click.option(
    "--submission-dir",
    metavar="PATH",
    type=click.Path(
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        writable=False,
        resolve_path=True,
    ),
    required=True,
    help="Path to the submission directory containing 'metadata/', 'files/', 'encrypted_files/' and 'logs/' directories",
)
@click.option(
    "--config-file",
    metavar="STRING",
    type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
    required=False,
    help="Path to config file",
)
@click.option(
    "--threads",
    default=min(len(sched_getaffinity(0)), 4),
    type=int,
    show_default=True,
    help="Number of threads to use for parallel operations",
)
def upload(
    submission_dir,
    config_file,
    threads,
):
    """
    Upload a submission to a GRZ/GDC.
    """
    from ..utils import read_config

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
