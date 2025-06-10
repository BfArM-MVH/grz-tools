"""Command for encrypting a submission."""

import logging
import sys
from pathlib import Path
from tempfile import NamedTemporaryFile

import click
from grz_common.cli import config_file, force, submission_dir
from grz_common.workers.worker import Worker

from ..models.config import EncryptConfig

log = logging.getLogger(__name__)


@click.command()
@submission_dir
@config_file
@force
def encrypt(submission_dir, config_file, force):
    """
    Encrypt a submission.

    Encryption is done with the recipient's public key.
    Sub-folders 'encrypted_files' and 'logs' are created within the submission directory.
    """
    config = EncryptConfig.from_path(config_file)

    submitter_privkey_path = config.keys.submitter_private_key_path
    if submitter_privkey_path == "":
        submitter_privkey_path = None

    log.info("Starting encryption...")

    submission_dir = Path(submission_dir)

    worker_inst = Worker(
        metadata_dir=submission_dir / "metadata",
        files_dir=submission_dir / "files",
        log_dir=submission_dir / "logs",
        encrypted_files_dir=submission_dir / "encrypted_files",
    )
    if pubkey := config.keys.grz_public_key:
        with NamedTemporaryFile("w") as f:
            f.write(pubkey)
            f.flush()
            worker_inst.encrypt(f.name, submitter_private_key_path=submitter_privkey_path, force=force)
    else:
        # This case cannot occur here, but an explicit check is needed for type-checking.
        if config.keys.grz_public_key_path is None:
            sys.exit("GRZ public key path is required for encryption.")
        worker_inst.encrypt(
            config.keys.grz_public_key_path,
            submitter_private_key_path=submitter_privkey_path,
            force=force,
        )

    log.info("Encryption successful!")
