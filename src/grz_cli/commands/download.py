"""Command for downloading a submission."""

import logging
from os import sched_getaffinity
from pathlib import Path

import click

from ..workers.worker import Worker

log = logging.getLogger(__name__)


@click.command()
@click.option(
    "--submission-id",
    required=True,
    type=str,
    metavar="STRING",
    help="S3 submission ID",
)
@click.option(
    "--output-dir",
    metavar="PATH",
    type=click.Path(
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        writable=True,
        resolve_path=True,
    ),
    required=True,
    default=None,
    help="Path to the target submission output directory",
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
def download(
    submission_id,
    output_dir,
    config_file,
    threads,
):
    """
    Download a submission from a GRZ.

    Downloaded metadata is stored within the `metadata` sub-folder of the submission output directory.
    Downloaded files are stored within the `encrypted_files` sub-folder of the submission output directory.
    """
    from ..utils import read_config

    config = read_config(config_file)

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
    worker_inst.download(config, submission_id)

    log.info("Download finished!")
