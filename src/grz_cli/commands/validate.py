"""Command for validating a submission."""

import logging
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
def validate(submission_dir):
    """
    Validate the submission.

    This validates the submission by checking its checksums, as well as performing basic sanity checks on the supplied metadata.
    Must be executed before calling `encrypt` and `upload`.
    """
    log.info("Starting validation...")

    submission_dir = Path(submission_dir)

    worker_inst = Worker(
        metadata_dir=submission_dir / "metadata",
        files_dir=submission_dir / "files",
        log_dir=submission_dir / "logs",
        encrypted_files_dir=submission_dir / "encrypted_files",
    )
    worker_inst.validate()

    log.info("Validation finished!")
