"""grzctl-specific wrappers for grz-cli commands that use the unified config."""

import logging
from pathlib import Path
from typing import TYPE_CHECKING

import click
import grz_cli.commands.encrypt as encrypt_module
import grz_cli.commands.validate as validate_module
import grz_common.cli as grzcli
from grz_common.workers.worker import Worker
from grz_db.errors import DuplicateInitialSubmissionError
from grz_db.models.submission import SubmissionStateEnum
from grz_pydantic_models.submission.metadata import GrzSubmissionMetadata

from ..commands import grzctl_configuration, inbox_option
from ..commands.db.cli import get_submission_db_instance
from ..dbcontext import DbContext
from ..models.config import GrzctlConfig

if TYPE_CHECKING:
    from grz_db.models.submission import SubmissionDb

log = logging.getLogger(__name__)


def derive_validate_config(config: GrzctlConfig) -> dict:
    """Build the config dict expected by grz-cli's ``validate`` command.

    ``ValidateConfig`` only needs ``identifiers`` (specifically ``identifiers.grz``).
    """
    return {"identifiers": config.identifiers.model_dump(mode="json", exclude_none=True)}


def derive_encrypt_config(config: GrzctlConfig) -> dict:
    """Build the config dict expected by grz-cli's ``encrypt`` command.

    ``EncryptConfig`` needs ``keys`` with ``grz_public_key_path`` (required,
    exactly one of ``grz_public_key`` or ``grz_public_key_path`` must be set)
    and optionally ``submitter_private_key_path`` and ``grz_private_key_path``.
    """
    keys_dump: dict = {}
    if config.keys.grz_public_key_path is not None:
        keys_dump["grz_public_key_path"] = config.keys.grz_public_key_path
    if config.keys.submitter_private_key_path is not None:
        keys_dump["submitter_private_key_path"] = config.keys.submitter_private_key_path
    if config.keys.grz_private_key_path is not None:
        keys_dump["grz_private_key_path"] = config.keys.grz_private_key_path
    return {"keys": keys_dump}


def _check_duplicate_initial(db: "SubmissionDb", metadata: GrzSubmissionMetadata) -> None:
    """Raise :class:`DuplicateInitialSubmissionError` if this initial submission is a duplicate.

    The rule and the index enforcing it belong to the database layer, which answers this in
    one transaction; all that is left here is unpacking the metadata. Validation is what the
    answer buys: the database would reject the submission anyway, but only once basic QC is
    being recorded, by which point the effort has been spent.

    A resolution failure is left to propagate. Reporting "not a duplicate" when the
    question could not be answered is how a duplicate would pass basic QC unnoticed, and
    nothing downstream treats an unresolvable case as an error any more.
    """
    submission = metadata.submission
    db.assert_no_duplicate_initial(
        metadata.submission_id,
        submitter_id=submission.submitter_id,
        local_case_id=submission.local_case_id,
        submission_type=submission.submission_type,
    )


def _warn_on_duplicate_initial(configuration: GrzctlConfig, submission_dir: Path, threads: int) -> None:
    """Best-effort, read-only duplicate-initial check when DB updates are disabled.

    Warns instead of raising, since nothing is written to the database in this mode.

    A check that was attempted and could not be answered is worth saying out loud:
    silence here reads as "not a duplicate", and that is an answer nobody got.
    """
    try:
        db = get_submission_db_instance(configuration.db.database_url, author=None)
        metadata = (
            Worker(
                metadata_dir=submission_dir / "metadata",
                files_dir=submission_dir / "files",
                log_dir=submission_dir / "logs",
                encrypted_files_dir=submission_dir / "encrypted_files",
                threads=threads,
            )
            .parse_submission()
            .metadata.content
        )
        _check_duplicate_initial(db, metadata)
    except DuplicateInitialSubmissionError as e:
        log.warning(f"{e} Basic QC would be recorded as failed; continuing because DB updates are disabled.")
    except Exception as e:
        log.warning(f"Could not check whether this is a duplicate initial submission: {e}")


@click.command()
@grzctl_configuration
@grzcli.submission_dir
@grzcli.force
@grzcli.threads
@click.option(
    "--mmap/--no-mmap",
    "mmap",
    default=False,
    hidden=True,
    help="Whether to use mmap.",
)
@grzcli.update_db
def validate(  # noqa: PLR0913
    configuration: GrzctlConfig,
    submission_dir,
    force,
    threads,
    mmap,
    update_db,
    **kwargs,
):
    """Validate the submission (wrapper with DB updates)."""
    submission_dir = Path(submission_dir)

    submission_id = ""
    if update_db:
        worker_inst = Worker(
            metadata_dir=submission_dir / "metadata",
            files_dir=submission_dir / "files",
            log_dir=submission_dir / "logs",
            encrypted_files_dir=submission_dir / "encrypted_files",
            threads=threads,
        )
        submission = worker_inst.parse_submission()
        submission_id = submission.metadata.content.submission_id

    with DbContext(
        configuration=configuration,
        submission_id=submission_id,
        start_state=SubmissionStateEnum.VALIDATING,
        end_state=SubmissionStateEnum.VALIDATED,
        enabled=update_db,
    ) as dbcontext_inst:
        if update_db:
            try:
                _check_duplicate_initial(dbcontext_inst.db, submission.metadata.content)
            except DuplicateInitialSubmissionError as e:
                # Another initial submission of this case already passed basic QC, so this
                # submission fails basic QC without spending any validation effort. Re-raising
                # lets DbContext record the ERROR state with the duplicate_initial reason.
                log.warning(f"{e} Failing basic QC for submission '{submission_id}' without validating.")
                _ = dbcontext_inst.db.modify_submission(submission_id, "basic_qc_passed", "false")
                raise
        else:
            _warn_on_duplicate_initial(configuration, submission_dir, threads)

        validate_module.validate.callback(  # type: ignore[misc]
            configuration=derive_validate_config(configuration),
            submission_dir=submission_dir,
            force=force,
            threads=threads,
            mmap=mmap,
            **kwargs,
        )

        if update_db:
            try:
                _ = dbcontext_inst.db.modify_submission(submission_id, "basic_qc_passed", "true")
            except DuplicateInitialSubmissionError as e:
                # Catches a competing initial submission of this case that passed basic QC
                # while this submission was being validated; handled the same way as the
                # pre-check above.
                log.warning(
                    f"Submission '{submission_id}' data validated, but {e} Failing basic QC for this submission."
                )
                _ = dbcontext_inst.db.modify_submission(submission_id, "basic_qc_passed", "false")
                raise


@click.command()
@grzctl_configuration
@grzcli.submission_dir
@grzcli.force
@click.option(
    "--check-validation-logs/--no-check-validation-logs",
    "check_validation_logs",
    default=True,
    help="Check validation logs before encrypting.",
)
@grzcli.update_db
def encrypt(
    configuration: GrzctlConfig,
    submission_dir,
    force,
    check_validation_logs,
    update_db,
    **kwargs,
):
    """Encrypt a submission (wrapper with DB updates)."""
    submission_dir = Path(submission_dir)

    submission_id = "unknown"
    if update_db:
        worker_inst = Worker(
            metadata_dir=submission_dir / "metadata",
            files_dir=submission_dir / "files",
            log_dir=submission_dir / "logs",
            encrypted_files_dir=submission_dir / "encrypted_files",
        )
        submission = worker_inst.parse_submission()
        submission_id = submission.metadata.content.submission_id

    with DbContext(
        configuration=configuration,
        submission_id=submission_id,
        start_state=SubmissionStateEnum.ENCRYPTING,
        end_state=SubmissionStateEnum.ENCRYPTED,
        enabled=update_db,
    ):
        encrypt_module.encrypt.callback(  # type: ignore[misc]
            configuration=derive_encrypt_config(configuration),
            submission_dir=submission_dir,
            force=force,
            check_validation_logs=check_validation_logs,
            **kwargs,
        )


@click.command()
@grzcli.submission_dir
@grzcli.threads
@grzctl_configuration
@inbox_option
@grzcli.update_db
def upload(
    configuration: GrzctlConfig,
    submission_dir,
    threads: int,
    inbox_name: str,
    update_db: bool,
    **kwargs,
):
    """Upload a submission to a GRZ/GDC (wrapper with DB updates)."""
    submission_dir = Path(submission_dir)

    worker_inst = Worker(
        metadata_dir=submission_dir / "metadata",
        files_dir=submission_dir / "files",
        log_dir=submission_dir / "logs",
        encrypted_files_dir=submission_dir / "encrypted_files",
        threads=threads,
    )
    submission = worker_inst.parse_submission()
    submission_id = submission.metadata.content.submission_id
    submitter_id = submission.metadata.content.submission.submitter_id

    s3_options = configuration.resolve_inbox(submitter_id=submitter_id, inbox_name=inbox_name).s3

    with DbContext(
        configuration=configuration,
        submission_id=submission_id,
        start_state=SubmissionStateEnum.UPLOADING,
        end_state=SubmissionStateEnum.UPLOADED,
        enabled=update_db,
    ):
        uploaded_id = worker_inst.upload(s3_options)
        click.echo(uploaded_id)
