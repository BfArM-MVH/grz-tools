"""Command for decrypting a submission."""

import logging
from pathlib import Path
from typing import Literal

import click
import grz_common.cli as grzcli
from grz_common.exceptions import ConfigurationError
from grz_common.utils.crypt import Crypt4GH
from grz_common.workers.worker import Worker
from grz_db.models.submission import SubmissionStateEnum

from ..commands import grzctl_configuration, inbox_option
from ..dbcontext import DbContext
from ..models.config import GrzctlConfig
from .inbox_resolution import require_inbox

log = logging.getLogger(__name__)

_ARCHIVES: dict[str, Literal["consented", "non_consented"]] = {
    "consented": "consented",
    "non-consented": "non_consented",
}
"""Maps each value of --archive to the name of the archive in the config."""


@click.command()
@grzctl_configuration
@grzcli.submission_dir
@inbox_option
@click.option(
    "--archive",
    type=click.Choice(list(_ARCHIVES)),
    default=None,
    help="Decrypt an archived submission with the private key of this archive from the config.",
)
@click.option(
    "--private-key-path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help=(
        "Decrypt with this crypt4gh private key instead of a key from the config. "
        "Its passphrase comes from C4GH_PASSPHRASE, else a prompt."
    ),
)
@grzcli.force
@grzcli.update_db
def decrypt(  # noqa: PLR0913, PLR0917
    configuration: GrzctlConfig,
    submission_dir,
    inbox_name,
    archive: str | None,
    private_key_path: Path | None,
    force,
    update_db,
    **kwargs,
):
    """
    Decrypt a submission.

    Decrypting a submission requires the _private_ key of the original recipient.
    By default, that is the private key of the inbox that the submission came from.
    grzctl looks up that inbox under the submitter named in the submission's metadata.
    The inbox is the one that --inbox names.
    Without --inbox, grzctl takes the inbox recorded in the database, else the submitter's only inbox.
    It reads the database only with --update-db.
    `grzctl download` records the inbox, and `grzctl db backfill` records it for older submissions.

    An archived submission is encrypted for its archive.
    --archive decrypts it with the private key of that archive, from archives.<archive>.private_key[_path].
    --private-key-path names a private key file, and that key wins over any key from the config.
    With --archive or --private-key-path, grzctl resolves no inbox, so neither goes with --inbox.

    --no-update-db keeps the state of the submission unchanged, for example when decrypting an archived submission.
    """
    if archive is not None and inbox_name is not None:
        raise click.UsageError("--archive and --inbox are mutually exclusive.")
    if private_key_path is not None and inbox_name is not None:
        raise click.UsageError("--private-key-path and --inbox are mutually exclusive.")

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
    ) as db_context:
        # raise a ConfigurationError, so that the DbContext records the failure reason configuration_error
        if private_key_path is not None:
            private_key = Crypt4GH.retrieve_private_key(private_key_path)
        elif archive is not None:
            private_key = configuration.archives.load_private_key(_ARCHIVES[archive])
        else:
            submitter_id = encrypted_submission.metadata.content.submission.submitter_id
            inbox_name = require_inbox(
                configuration,
                submitter_id=submitter_id,
                submission_id=submission_id,
                inbox_name=inbox_name,
                db_service=db_context.db,
                hint="Pass --inbox, or record the inbox with 'grzctl db backfill'.",
                exc_type=ConfigurationError,
            )
            private_key = configuration.inbox_target(submitter_id, inbox_name).load_private_key()
        worker_inst.decrypt(recipient_private_key=private_key, force=force)

    log.info("Decryption successful!")
