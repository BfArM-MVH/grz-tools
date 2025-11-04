"""Command for validating a submission."""

import logging
from pathlib import Path

import click
from grz_cli.models.config import ValidateConfig
from grz_common.cli import (
    config_file,
    files_dir,
    force,
    logs_dir,
    metadata_dir,
    submission_dir,
    threads,
)
from grz_common.workers.worker import Worker

log = logging.getLogger(__name__)


@click.command()
@submission_dir
@metadata_dir
@files_dir
@logs_dir
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
def validate(submission_dir, metadata_dir, files_dir, logs_dir, config_file, force, threads, with_grz_check):
    """
    Validate the submission.

    This validates the submission by checking its checksums, as well as performing basic sanity checks on the supplied metadata.
    Must be executed before calling `encrypt` and `upload`.
    """
    bundled_mode = submission_dir is not None
    granular_mode = any(map(lambda v: v is not None, [metadata_dir, files_dir, logs_dir]))

    if bundled_mode and granular_mode:
        raise click.UsageError("'--submission-dir' is mutually exclusive with explicit path options.")

    if bundled_mode:
        base = Path(submission_dir)
        _metadata_dir = base / "metadata"
        _files_dir = base / "files"
        _logs_dir = base / "logs"
    elif granular_mode:
        required = {
            "--metadata-dir": metadata_dir,
            "--files-dir": files_dir,
            "--logs-dir": logs_dir,
        }
        missing = [name for name, path in required.items() if path is None]
        if missing:
            raise click.UsageError(f"Flexible mode requires: {', '.join(missing)}")
        _metadata_dir, _files_dir, _logs_dir = Path(metadata_dir), Path(files_dir), Path(logs_dir)
    else:
        raise click.UsageError("You must specify either '--submission-dir' or the required explicit path options.")

    config = ValidateConfig.from_path(config_file)

    log.info("Starting validation...")

    _logs_dir.mkdir(parents=True, exist_ok=True)

    worker_inst = Worker(
        metadata_dir=_metadata_dir,
        files_dir=_files_dir,
        log_dir=_logs_dir,
        # encrypted_files_dir is not used by validate, but required by Worker
        encrypted_files_dir=_files_dir.parent / "encrypted_files",
        threads=threads,
    )
    worker_inst.validate(identifiers=config.identifiers, force=force, with_grz_check=with_grz_check)

    log.info("Validation finished!")
