"""Command for validating a submission."""

import logging
from pathlib import Path

import click
from grz_cli.models.config import ValidateConfig
from grz_common.cli import config_file, force, submission_dir, threads
from grz_common.workers.worker import Worker

log = logging.getLogger(__name__)


@click.command()
@submission_dir
@config_file
@force
@threads
@click.option(
    "--with-grz-check/--no-grz-check",
    "with_grz_check",
    default=True,
    hidden=True,
    help="Whether to use grz-check to perform validation",
)
def validate(submission_dir, config_file, force, threads, with_grz_check):
    """
    Validate the submission.

    This validates the submission by checking its checksums, as well as performing basic sanity checks on the supplied metadata.
    Must be executed before calling `encrypt` and `upload`.
    """
    config = ValidateConfig.from_path(config_file)

    log.info("Starting validation...")

    submission_dir = Path(submission_dir)

    worker_inst = Worker(
        metadata_dir=submission_dir / "metadata",
        files_dir=submission_dir / "files",
        log_dir=submission_dir / "logs",
        encrypted_files_dir=submission_dir / "encrypted_files",
        threads=threads,
    )
    worker_inst.validate(identifiers=config.identifiers, force=force, with_grz_check=with_grz_check)

    log.info("Validation finished!")
