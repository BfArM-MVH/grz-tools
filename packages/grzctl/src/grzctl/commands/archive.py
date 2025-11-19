"""Command for archiving a submission."""

import logging
from pathlib import Path

import click
from grz_common.cli import DIR_RW_C, config_file, logs_dir, submission_dir, threads
from grz_common.workers.worker import Worker

from ..models.config import ArchiveConfig

log = logging.getLogger(__name__)


@click.command()
@submission_dir
@click.option("--metadata-dir", type=DIR_RW_C, required=False)
@click.option("--encrypted-files-dir", type=DIR_RW_C, required=False)
@logs_dir
@config_file
@threads
def archive(submission_dir, metadata_dir, encrypted_files_dir, logs_dir, config_file, threads):  # noqa: PLR0913
    """
    Archive a pre-staged submission directory from within a GRZ/GDC.
    This command expects a directory containing redacted metadata, redacted logs,
    and encrypted files.
    """
    bundled_mode = submission_dir is not None
    granular_mode = any([metadata_dir, encrypted_files_dir, logs_dir])

    if bundled_mode and granular_mode:
        raise click.UsageError("'--output-dir' is mutually exclusive with explicit path options.")

    if bundled_mode:
        base = Path(submission_dir)
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
            raise click.UsageError(f"Flexible mode requires: {', '.join(missing)}")
        _metadata_dir, _encrypted_files_dir, _logs_dir = Path(metadata_dir), Path(encrypted_files_dir), Path(logs_dir)
    else:
        raise click.UsageError("You must specify either '--submission-dir' or all explicit path options.")

    config = ArchiveConfig.from_path(config_file)
    log.info("Starting archival upload...")

    worker_inst = Worker(
        metadata_dir=_metadata_dir,
        files_dir="/dev/null",
        log_dir=_logs_dir,
        encrypted_files_dir=_encrypted_files_dir,
        threads=threads,
    )
    worker_inst.archive(config.s3)

    log.info("Archival finished!")
