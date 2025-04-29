"""Command for encrypting a submission."""

import logging
from pathlib import Path
from tempfile import NamedTemporaryFile

import click

from ..utils.config import read_config
from ..workers.worker import Worker
from .common import config_file, submission_dir

log = logging.getLogger(__name__)


@click.command()
@submission_dir
@config_file
def encrypt(
    submission_dir,
    config_file,
):
    """
    Encrypt a submission.

    Encryption is done with the recipient's public key.
    Sub-folders 'encrypted_files' and 'logs' are created within the submission directory.
    """
    config = read_config(config_file)

    submitter_privkey_path = config.submitter_private_key_path
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
    if config.grz_public_key:
        with NamedTemporaryFile("w") as f:
            f.write(config.grz_public_key)
            f.flush()
            worker_inst.encrypt(
                f.name,
                submitter_private_key_path=submitter_privkey_path,
            )
    else:
        worker_inst.encrypt(
            config.grz_public_key_path,
            submitter_private_key_path=submitter_privkey_path,
        )

    log.info("Encryption successful!")
