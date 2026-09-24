"""Command for decrypting a submission."""

import logging
from collections.abc import Iterable
from pathlib import Path

import click
import grz_common.cli as grzcli
from grz_common.utils.crypt import Crypt4GH
from grz_common.workers.worker import Worker
from grz_db.models.submission import SubmissionStateEnum

from ..commands import grzctl_configuration
from ..dbcontext import DbContext
from ..models.config import GrzctlConfig

log = logging.getLogger(__name__)


@click.command()
@grzctl_configuration
@grzcli.submission_dir
@click.option(
    "--private-key-file",
    "private_key_file",
    type=click.Path(exists=True, dir_okay=False),
    default=None,
    help="Decrypt with this crypt4gh private key instead of the keys in the config. "
    "Its passphrase comes from C4GH_PASSPHRASE, else a prompt.",
)
@grzcli.force
@grzcli.update_db
def decrypt(
    configuration: GrzctlConfig,
    submission_dir,
    private_key_file,
    force,
    update_db,
    **kwargs,
):
    """
    Decrypt a submission.

    Decrypting a submission requires the _private_ key of the original recipient.
    Without --private-key-file, the keys in the config are tried one at a time:
    first the private keys of all inboxes of the submitter named in the submission's metadata,
    then the private keys of the consented and the non-consented archive.
    The keys are tested against the crypt4gh header of the submission's first encrypted file.
    The first key that opens it decrypts all files.
    """
    log.info("Starting decryption...")

    submission_dir = Path(submission_dir)

    worker_inst = Worker(
        metadata_dir=submission_dir / "metadata",
        files_dir=submission_dir / "files",
        log_dir=submission_dir / "logs",
        encrypted_files_dir=submission_dir / "encrypted_files",
    )
    encrypted_submission = worker_inst.parse_encrypted_submission()
    submission_id = encrypted_submission.submission_id

    with DbContext(
        configuration=configuration,
        submission_id=submission_id,
        start_state=SubmissionStateEnum.DECRYPTING,
        end_state=SubmissionStateEnum.DECRYPTED,
        enabled=update_db,
    ):
        private_keys: Iterable[tuple[str, bytes]]
        if private_key_file is not None:
            private_keys = [(f"--private-key-file {private_key_file}", Crypt4GH.retrieve_private_key(private_key_file))]
        else:
            submitter_id = encrypted_submission.metadata.content.submission.submitter_id
            private_keys = configuration.iter_decryption_keys(submitter_id)
        private_key = encrypted_submission.find_private_key(private_keys)
        worker_inst.decrypt(recipient_private_key=private_key, force=force)

    log.info("Decryption successful!")
