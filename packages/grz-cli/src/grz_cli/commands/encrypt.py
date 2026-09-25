"""Command for encrypting a submission."""

import logging
import sys
from pathlib import Path
from typing import Any

import click
import grz_common.cli as grzcli
from grz_common.utils.crypt import Crypt4GH
from grz_common.workers.worker import Worker

from ..models.config import EncryptConfig

log = logging.getLogger(__name__)


@click.command()
@grzcli.configuration
@grzcli.submission_dir
@grzcli.force
@click.option(
    "--check-validation-logs/--no-check-validation-logs",
    "check_validation_logs",
    default=True,
    help="Check validation logs before encrypting.",
)
def encrypt(configuration: dict[str, Any], submission_dir, force, check_validation_logs, **kwargs):
    """
    Encrypt a submission.

    Encryption is done with the recipient's public key.
    Sub-folders 'encrypted_files' and 'logs' are created within the submission directory.
    """
    config = EncryptConfig.model_validate(configuration)

    if config.keys.grz_public_key is not None:
        grz_public_key = Crypt4GH.load_public_key(config.keys.grz_public_key, key_name="keys.grz_public_key")
    elif config.keys.grz_public_key_path is not None:
        grz_public_key = Crypt4GH.retrieve_public_key(config.keys.grz_public_key_path)
    else:
        # This case cannot occur here, but an explicit check is needed for type-checking.
        sys.exit("Either keys.grz_public_key or keys.grz_public_key_path must be set for encryption.")

    submitter_private_key = None
    if config.keys.submitter_private_key is not None:
        submitter_private_key = Crypt4GH.load_private_key(
            config.keys.submitter_private_key.get_secret_value(), key_name="keys.submitter_private_key"
        )
    elif config.keys.submitter_private_key_path is not None:
        submitter_private_key = Crypt4GH.retrieve_private_key(config.keys.submitter_private_key_path)

    log.info("Starting encryption...")

    submission_dir = Path(submission_dir)

    worker_inst = Worker(
        metadata_dir=submission_dir / "metadata",
        files_dir=submission_dir / "files",
        log_dir=submission_dir / "logs",
        encrypted_files_dir=submission_dir / "encrypted_files",
    )
    worker_inst.encrypt(
        grz_public_key,
        submitter_private_key=submitter_private_key,
        force=force,
        check_validation_logs=check_validation_logs,
    )

    log.info("Encryption successful!")
