"""Command for encrypting a submission."""

import logging
import sys
from pathlib import Path
from typing import Any

import click
import grz_common.cli as grzcli
from grz_common.models.base import get_secret_value
from grz_common.utils.crypt import Crypt4GH
from grz_common.workers.worker import Worker

from ..models.config import EncryptConfig

log = logging.getLogger(__name__)


@click.command()
@grzcli.configuration
@grzcli.submission_dir
@grzcli.metadata_dir
@grzcli.files_dir
@grzcli.output_encrypted_files_dir
@grzcli.logs_dir
@grzcli.force
@click.option(
    "--check-validation-logs/--no-check-validation-logs",
    "check_validation_logs",
    default=True,
    help="Check validation logs before encrypting.",
)
def encrypt(
    configuration: dict[str, Any],
    submission_dir,
    metadata_dir,
    files_dir,
    output_encrypted_files_dir,
    logs_dir,
    force,
    check_validation_logs,
    **kwargs,
):
    """
    Encrypt a submission.

    Encryption is done with the recipient's public key.
    """
    bundled_mode = submission_dir is not None
    granular_mode = any(map(lambda v: v is not None, [metadata_dir, files_dir, output_encrypted_files_dir, logs_dir]))

    if bundled_mode and granular_mode:
        raise click.UsageError("'--submission-dir' is mutually exclusive with explicit path options.")

    if bundled_mode:
        base = Path(submission_dir)
        _metadata_dir = base / "metadata"
        _files_dir = base / "files"
        _encrypted_files_dir = base / "encrypted_files"
        _logs_dir = base / "logs"
    elif granular_mode:
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

    config = EncryptConfig.model_validate(configuration)

    if config.keys.grz_public_key is not None:
        grz_public_key = Crypt4GH.load_public_key(config.keys.grz_public_key, key_name="keys.grz_public_key")
    elif config.keys.grz_public_key_path is not None:
        grz_public_key = Crypt4GH.retrieve_public_key(config.keys.grz_public_key_path)
    else:
        # This case cannot occur here, but an explicit check is needed for type-checking.
        sys.exit("Either keys.grz_public_key or keys.grz_public_key_path must be set for encryption.")

    submitter_private_key = None
    passphrase = get_secret_value(config.keys.submitter_private_key_passphrase)
    if config.keys.submitter_private_key is not None:
        submitter_private_key = Crypt4GH.load_private_key(
            config.keys.submitter_private_key.get_secret_value(),
            passphrase=passphrase,
            key_name="keys.submitter_private_key",
        )
    elif config.keys.submitter_private_key_path is not None:
        submitter_private_key = Crypt4GH.retrieve_private_key(
            config.keys.submitter_private_key_path, passphrase=passphrase
        )

    log.info("Starting encryption...")

    _encrypted_files_dir.parent.mkdir(parents=True, exist_ok=True)
    _logs_dir.parent.mkdir(parents=True, exist_ok=True)

    worker_inst = Worker(
        metadata_dir=_metadata_dir,
        files_dir=_files_dir,
        log_dir=_logs_dir,
        encrypted_files_dir=_encrypted_files_dir,
    )
    worker_inst.encrypt(
        grz_public_key,
        submitter_private_key=submitter_private_key,
        force=force,
        check_validation_logs=check_validation_logs,
    )

    log.info("Encryption successful!")
