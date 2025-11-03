"""Command for decrypting a submission."""

import logging
import sys
from pathlib import Path

import click
from grz_common.cli import (
    config_file,
    encrypted_files_dir,
    force,
    metadata_dir,
    output_files_dir,
    output_logs_dir,
    submission_dir,
)
from grz_common.workers.worker import Worker

from ..models.config import DecryptConfig

log = logging.getLogger(__name__)


@click.command()
@submission_dir
@metadata_dir
@encrypted_files_dir
@output_files_dir
@output_logs_dir
@config_file
@force
def decrypt(submission_dir, metadata_dir, encrypted_files_dir, output_files_dir, output_logs_dir, config_file, force):
    """
    Decrypt a submission.

    Decrypting a submission requires the _private_ key of the original recipient.
    """
    in_legacy_mode = submission_dir is not None
    in_flexible_mode = any([metadata_dir, encrypted_files_dir, output_files_dir, output_logs_dir])

    if in_legacy_mode and in_flexible_mode:
        raise click.UsageError("'--submission-dir' is mutually exclusive with explicit path options.")

    if in_legacy_mode:
        log.warning("Using deprecated '--submission-dir'. Please switch to explicit path options.")
        base = Path(submission_dir)
        _metadata_dir = base / "metadata"
        _encrypted_files_dir = base / "encrypted_files"
        _files_dir = base / "files"
        _logs_dir = base / "logs"
    elif in_flexible_mode:
        required = {
            "--metadata-dir": metadata_dir,
            "--encrypted-files-dir": encrypted_files_dir,
            "--output-files-dir": output_files_dir,
            "--output-logs-dir": output_logs_dir,
        }
        missing = [name for name, path in required.items() if path is None]
        if missing:
            raise click.UsageError(f"Flexible mode requires: {', '.join(missing)}")
        _metadata_dir, _encrypted_files_dir, _files_dir, _logs_dir = Path(metadata_dir), Path(encrypted_files_dir), Path(output_files_dir), Path(output_logs_dir)
    else:
        raise click.UsageError("You must specify either '--submission-dir' or the required explicit path options.")

    config = DecryptConfig.from_path(config_file)

    grz_privkey_path = config.keys.grz_private_key_path
    if not grz_privkey_path:
        log.error("GRZ private key path is required for decryption.")
        sys.exit(1)

    log.info("Starting decryption...")

    _files_dir.parent.mkdir(parents=True, exist_ok=True)
    _logs_dir.parent.mkdir(parents=True, exist_ok=True)

    worker_inst = Worker(
        metadata_dir=_metadata_dir,
        files_dir=_files_dir,
        log_dir=_logs_dir,
        encrypted_files_dir=_encrypted_files_dir,
    )
    worker_inst.decrypt(grz_privkey_path, force=force)

    log.info("Decryption successful!")
