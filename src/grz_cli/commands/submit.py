"""Command for submitting (validating, encrypting, and uploading) a submission."""

import logging
from os import sched_getaffinity

import click

log = logging.getLogger(__name__)


@click.command("submit")
@click.option(
    "--submission-dir",
    metavar="PATH",
    type=click.Path(
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        writable=False,
        resolve_path=True,
    ),
    required=True,
    help="Path to the submission directory containing 'metadata/', 'files/', 'encrypted_files/' and 'logs/' directories",
)
@click.option(
    "--config-file",
    metavar="STRING",
    type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
    required=False,
    help="Path to config file",
)
@click.option(
    "--threads",
    default=min(len(sched_getaffinity(0)), 4),
    type=int,
    show_default=True,
    help="Number of threads to use for parallel operations",
)
@click.pass_context
def submit(ctx, submission_dir, config_file, threads):
    """
    Validate, encrypt, and then upload.

    This is a convenience command that performs the following steps in order:
    1. Validate the submission
    2. Encrypt the submission
    3. Upload the encrypted submission
    """
    from .encrypt import encrypt
    from .upload import upload
    from .validate import validate

    click.echo("Starting submission process...")
    ctx.invoke(validate, submission_dir=submission_dir)
    ctx.invoke(encrypt, submission_dir=submission_dir, config_file=config_file)
    ctx.invoke(upload, submission_dir=submission_dir, config_file=config_file, threads=threads)
    click.echo("Submission finished!")
