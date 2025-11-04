"""Command for encrypting a submission."""

import logging
import sys
from pathlib import Path
from tempfile import NamedTemporaryFile

import click
from grz_common.cli import (
    config_file,
    files_dir,
    force,
    logs_dir,
    metadata_dir,
    output_encrypted_files_dir,
    submission_dir,
)
from grz_common.workers.worker import Worker

from ..models.config import EncryptConfig

log = logging.getLogger(__name__)


@click.command()
@submission_dir
@metadata_dir
@files_dir
@output_encrypted_files_dir
@logs_dir
@config_file
@force
@click.option(
    "--check-validation-logs/--no-check-validation-logs",
    "check_validation_logs",
    default=True,
    help="Check validation logs before encrypting.",
)
def encrypt(
    submission_dir,
    metadata_dir,
    files_dir,
    output_encrypted_files_dir,
    logs_dir,
    config_file,
    force,
    check_validation_logs,
):
    """
    Encrypt a submission.

    Encryption is done with the recipient's public key.
    """
    in_legacy_mode = submission_dir is not None
    in_flexible_mode = any(
        map(lambda v: v is not None, [metadata_dir, files_dir, output_encrypted_files_dir, logs_dir])
    )

    if in_legacy_mode and in_flexible_mode:
        raise click.UsageError("'--submission-dir' is mutually exclusive with explicit path options.")

    if in_legacy_mode:
        base = Path(submission_dir)
        _metadata_dir = base / "metadata"
        _files_dir = base / "files"
        _encrypted_files_dir = base / "encrypted_files"
        _logs_dir = base / "logs"
    elif in_flexible_mode:
        required = {
            "--metadata-dir": metadata_dir,
            "--files-dir": files_dir,
            "--logs-dir": logs_dir,
            "--output-encrypted-files-dir": output_encrypted_files_dir,
        }
        missing = [name for name, path in required.items() if path is None]
        if missing:
            raise click.UsageError(f"Flexible mode requires: {', '.join(missing)}")
        _metadata_dir, _files_dir, _encrypted_files_dir, _logs_dir = (
            Path(metadata_dir),
            Path(files_dir),
            Path(output_encrypted_files_dir),
            Path(logs_dir),
        )
    else:
        raise click.UsageError("You must specify either '--submission-dir' or the required explicit path options.")

    config = EncryptConfig.from_path(config_file)

    submitter_privkey_path = config.keys.submitter_private_key_path
    if submitter_privkey_path == "":
        submitter_privkey_path = None

    log.info("Starting encryption...")

    _encrypted_files_dir.parent.mkdir(parents=True, exist_ok=True)
    _logs_dir.parent.mkdir(parents=True, exist_ok=True)

    worker_inst = Worker(
        metadata_dir=_metadata_dir,
        files_dir=_files_dir,
        log_dir=_logs_dir,
        encrypted_files_dir=_encrypted_files_dir,
    )
    if pubkey := config.keys.grz_public_key:
        with NamedTemporaryFile("w") as f:
            f.write(pubkey)
            f.flush()
            worker_inst.encrypt(
                f.name,
                submitter_private_key_path=submitter_privkey_path,
                force=force,
                check_validation_logs=check_validation_logs,
            )
    else:
        # This case cannot occur here, but an explicit check is needed for type-checking.
        if config.keys.grz_public_key_path is None:
            sys.exit("GRZ public key path is required for encryption.")
        worker_inst.encrypt(
            config.keys.grz_public_key_path,
            submitter_private_key_path=submitter_privkey_path,
            force=force,
            check_validation_logs=check_validation_logs,
        )

    log.info("Encryption successful!")
