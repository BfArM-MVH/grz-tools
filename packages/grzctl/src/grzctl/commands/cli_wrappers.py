"""grzctl-specific wrappers for grz-cli commands that use the unified config."""

import logging
from pathlib import Path

import click
import grz_cli.commands.encrypt as encrypt_module
import grz_cli.commands.validate as validate_module
import grz_common.cli as grzcli
from grz_common.workers.worker import Worker
from grz_db.models.submission import SubmissionStateEnum

from ..commands import grzctl_configuration, inbox_option
from ..dbcontext import DbContext
from ..models.config import GrzctlConfig

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


@click.command()
@grzctl_configuration
@grzcli.submission_dir
@grzcli.metadata_dir
@grzcli.files_dir
@grzcli.logs_dir
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
    metadata_dir,
    files_dir,
    logs_dir,
    force,
    threads,
    mmap,
    update_db,
    **kwargs,
):
    """Validate the submission (wrapper with DB updates)."""
    bundled_mode = submission_dir is not None
    granular_mode = any(map(lambda v: v is not None, [metadata_dir, files_dir, logs_dir]))

    if bundled_mode and granular_mode:
        raise click.UsageError("'--submission-dir' is mutually exclusive with explicit path options.")

    if not bundled_mode and not granular_mode:
        raise click.UsageError("You must specify either '--submission-dir' or the required explicit path options.")

    _metadata_dir = Path(submission_dir) / "metadata" if bundled_mode else Path(metadata_dir)
    _files_dir = Path(submission_dir) / "files" if bundled_mode else Path(files_dir)
    _logs_dir = Path(submission_dir) / "logs" if bundled_mode else Path(logs_dir)
    _encrypted_files_dir = (
        Path(submission_dir) / "encrypted_files" if bundled_mode else Path(files_dir).parent / "encrypted_files"
    )

    worker_inst = Worker(
        metadata_dir=_metadata_dir,
        files_dir=_files_dir,
        log_dir=_logs_dir,
        encrypted_files_dir=_encrypted_files_dir,
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
        validate_module.validate.callback(  # type: ignore[misc]
            configuration=derive_validate_config(configuration),
            submission_dir=submission_dir,
            metadata_dir=metadata_dir,
            files_dir=files_dir,
            logs_dir=logs_dir,
            force=force,
            threads=threads,
            mmap=mmap,
        )

        if update_db:
            _ = dbcontext_inst.db.modify_submission(submission_id, "basic_qc_passed", "true")


@click.command()
@grzctl_configuration
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
@grzcli.update_db
def encrypt(  # noqa: PLR0913
    configuration: GrzctlConfig,
    submission_dir,
    metadata_dir,
    files_dir,
    output_encrypted_files_dir,
    logs_dir,
    force,
    check_validation_logs,
    update_db,
    **kwargs,
):
    """Encrypt a submission (wrapper with DB updates)."""
    bundled_mode = submission_dir is not None
    granular_mode = any(map(lambda v: v is not None, [metadata_dir, files_dir, output_encrypted_files_dir, logs_dir]))

    if bundled_mode and granular_mode:
        raise click.UsageError("'--submission-dir' is mutually exclusive with explicit path options.")

    if not bundled_mode and not granular_mode:
        raise click.UsageError("You must specify either '--submission-dir' or the required explicit path options.")

    _metadata_dir = Path(submission_dir) / "metadata" if bundled_mode else Path(metadata_dir)
    _files_dir = Path(submission_dir) / "files" if bundled_mode else Path(files_dir)
    _logs_dir = Path(submission_dir) / "logs" if bundled_mode else Path(logs_dir)
    _encrypted_files_dir = (
        Path(submission_dir) / "encrypted_files" if bundled_mode else Path(output_encrypted_files_dir)
    )

    worker_inst = Worker(
        metadata_dir=_metadata_dir,
        files_dir=_files_dir,
        log_dir=_logs_dir,
        encrypted_files_dir=_encrypted_files_dir,
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
            submission_dir=submission_dir if bundled_mode else None,
            metadata_dir=metadata_dir,
            files_dir=files_dir,
            output_encrypted_files_dir=output_encrypted_files_dir,
            logs_dir=logs_dir,
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
