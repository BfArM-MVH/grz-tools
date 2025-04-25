"""
CLI module for handling command-line interface operations.
"""

import importlib
import json
import logging
import logging.config
import os
import sys
from os import PathLike, sched_getaffinity
from pathlib import Path
from tempfile import NamedTemporaryFile

import click
import grz_pydantic_models.submission.metadata
import platformdirs
import rich.console
import rich.table
import yaml
from pydantic_core import to_jsonable_python

from .constants import PACKAGE_ROOT
from .download import query_submissions
from .logging_setup import add_filelogger
from .models.config import ConfigModel
from .parser import Worker

log = logging.getLogger(PACKAGE_ROOT + ".cli")

DEFAULT_CONFIG_PATH = Path(platformdirs.user_config_dir("grz-cli")) / "config.yaml"

# Aliases for path types for click options
# Naming convention: {DIR,FILE}_{Read,Write}_{Exists,Create}
DIR_R_E = click.Path(
    exists=True,
    file_okay=False,
    dir_okay=True,
    readable=True,
    writable=False,
    resolve_path=True,
)
DIR_RW_E = click.Path(
    exists=True,
    file_okay=False,
    dir_okay=True,
    readable=True,
    writable=True,
    resolve_path=True,
)
DIR_RW_C = click.Path(
    exists=False,
    file_okay=False,
    dir_okay=True,
    readable=True,
    writable=True,
    resolve_path=True,
)
FILE_R_E = click.Path(exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True)

submission_dir = click.option(
    "--submission-dir",
    metavar="PATH",
    type=DIR_R_E,
    required=True,
    help="Path to the submission directory containing 'metadata/', 'files/', 'encrypted_files/' and 'logs/' directories",
)

config_file = click.option(
    "--config-file",
    metavar="STRING",
    type=FILE_R_E,
    required=False,
    default=DEFAULT_CONFIG_PATH,
    help="Path to config file",
)

threads = click.option(
    "--threads",
    default=min(len(sched_getaffinity(0)), 4),
    type=int,
    show_default=True,
    help="Number of threads to use for parallel operations",
)

submission_id = click.option(
    "--submission-id",
    required=True,
    type=str,
    metavar="STRING",
    help="S3 submission ID",
)

output_dir = click.option(
    "--output-dir",
    metavar="PATH",
    type=DIR_RW_E,
    required=True,
    default=None,
    help="Path to the target submission output directory",
)

output_json = click.option("--json", "output_json", is_flag=True, help="Output JSON for machine-readability.")


def read_config(config_path: str | PathLike) -> ConfigModel:
    """Reads the configuration file and validates it against the schema."""
    with open(config_path, encoding="utf-8") as f:
        config = ConfigModel(**yaml.safe_load(f))

    return config


class OrderedGroup(click.Group):
    """
    A click Group that keeps track of the order in which commands are added.
    """

    def list_commands(self, ctx):
        """Return the list of commands in the order they were added."""
        return list(self.commands.keys())


@click.command()
@submission_dir
def validate(submission_dir):
    """
    Validate the submission.

    This validates the submission by checking its checksums, as well as performing basic sanity checks on the supplied metadata.
    Must be executed before calling `encrypt` and `upload`.
    """
    log.info("Starting validation...")

    submission_dir = Path(submission_dir)

    worker_inst = Worker(
        metadata_dir=submission_dir / "metadata",
        files_dir=submission_dir / "files",
        log_dir=submission_dir / "logs",
        encrypted_files_dir=submission_dir / "encrypted_files",
    )
    worker_inst.validate()

    log.info("Validation finished!")


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


@click.command()
@submission_dir
@config_file
@threads
def upload(
    submission_dir,
    config_file,
    threads,
):
    """
    Upload a submission to a GRZ/GDC.
    """
    config = read_config(config_file)

    log.info("Starting upload...")

    submission_dir = Path(submission_dir)

    worker_inst = Worker(
        metadata_dir=submission_dir / "metadata",
        files_dir=submission_dir / "files",
        log_dir=submission_dir / "logs",
        encrypted_files_dir=submission_dir / "encrypted_files",
        threads=threads,
    )
    # output the generated submission ID
    print(worker_inst.upload(config))

    log.info("Upload finished!")


@click.command("submit")
@submission_dir
@config_file
@threads
@click.pass_context
def submit(ctx, submission_dir, config_file, threads):
    """
    Validate, encrypt, and then upload.

    This is a convenience command that performs the following steps in order:
    1. Validate the submission
    2. Encrypt the submission
    3. Upload the encrypted submission
    """
    click.echo("Starting submission process...")
    ctx.invoke(validate, submission_dir=submission_dir)
    ctx.invoke(encrypt, submission_dir=submission_dir, config_file=config_file)
    ctx.invoke(upload, submission_dir=submission_dir, config_file=config_file, threads=threads)
    click.echo("Submission finished!")


@click.command()
@submission_id
@output_dir
@config_file
@threads
def download(
    submission_id,
    output_dir,
    config_file,
    threads,
):
    """
    Download a submission from a GRZ.

    Downloaded metadata is stored within the `metadata` sub-folder of the submission output directory.
    Downloaded files are stored within the `encrypted_files` sub-folder of the submission output directory.
    """
    config = read_config(config_file)

    log.info("Starting download...")

    submission_dir_path = Path(output_dir)
    if not submission_dir_path.is_dir():
        log.debug("Creating submission directory %s", submission_dir_path)
        submission_dir_path.mkdir(mode=0o770, parents=False, exist_ok=False)

    worker_inst = Worker(
        metadata_dir=submission_dir_path / "metadata",
        files_dir=submission_dir_path / "files",
        log_dir=submission_dir_path / "logs",
        encrypted_files_dir=submission_dir_path / "encrypted_files",
        threads=threads,
    )
    worker_inst.download(config, submission_id)

    log.info("Download finished!")


@click.command()
@config_file
@output_json
def list_submissions(config_file, output_json):
    """
    List submissions within an inbox from oldest to newest.
    """
    config = read_config(config_file)
    submissions = query_submissions(config)

    if output_json:
        json.dump(to_jsonable_python(submissions), sys.stdout)
    else:
        console = rich.console.Console()
        table = rich.table.Table()
        table.add_column("ID", no_wrap=True)
        table.add_column("Status", no_wrap=True)
        table.add_column("Oldest Upload", overflow="fold")
        table.add_column("Newest Upload", overflow="fold")
        for submission in submissions:
            status_text = rich.text.Text(
                "Complete" if submission.complete else "Incomplete",
                style="green" if submission.complete else "yellow",
            )
            table.add_row(
                submission.submission_id,
                status_text,
                submission.oldest_upload.astimezone().strftime("%Y-%m-%d %H:%M:%S"),
                submission.newest_upload.astimezone().strftime("%Y-%m-%d %H:%M:%S"),
            )
        console.print(table)


@click.command()
@submission_dir
@config_file
def decrypt(
    submission_dir,
    config_file,
):
    """
    Decrypt a submission.

    Decrypting a submission requires the _private_ key of the original recipient.
    """
    config = read_config(config_file)

    grz_privkey_path = config.grz_private_key_path
    if not grz_privkey_path:
        log.error("GRZ private key path is required for decryption.")
        sys.exit(1)

    log.info("Starting decryption...")

    submission_dir = Path(submission_dir)

    worker_inst = Worker(
        metadata_dir=submission_dir / "metadata",
        files_dir=submission_dir / "files",
        log_dir=submission_dir / "logs",
        encrypted_files_dir=submission_dir / "encrypted_files",
    )
    worker_inst.decrypt(grz_privkey_path)

    log.info("Encryption successful!")


def build_cli(grz_mode=False):
    """
    Factory for building the CLI application.
    """

    @click.group(
        cls=OrderedGroup,
        help="Validate, encrypt, decrypt and upload submissions to a GRZ/GDC.",
    )
    @click.version_option(
        version=importlib.metadata.version("grz-cli"),
        prog_name="grz-cli",
        message=f"%(prog)s v%(version)s (metadata schema versions: {', '.join(grz_pydantic_models.submission.metadata.get_supported_versions())})",
    )
    @click.option("--log-file", metavar="FILE", type=str, help="Path to log file")
    @click.option(
        "--log-level",
        default="INFO",
        type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]),
        help="Set the log level (default: INFO)",
    )
    def cli(log_file: str | None = None, log_level: str = "INFO"):
        """
        Command-line interface function for setting up logging.

        :param log_file: Path to the log file. If provided, a file logger will be added.
        :param log_level: Log level for the logger. It should be one of the following:
                           DEBUG, INFO, WARNING, ERROR, CRITICAL.
        """
        if log_file:
            add_filelogger(
                log_file,
                log_level.upper(),
            )  # Add file logger

        # show only time and log level in STDOUT
        logging.basicConfig(
            format="%(asctime)s - %(levelname)s - %(message)s",
        )

        # set the log level for this package
        logging.getLogger(PACKAGE_ROOT).setLevel(log_level.upper())

        log.debug("Logging setup complete.")

    cli.add_command(validate)
    cli.add_command(encrypt)
    cli.add_command(upload)
    cli.add_command(submit)

    if grz_mode:
        cli.add_command(list_submissions, name="list")
        cli.add_command(download)
        cli.add_command(decrypt)

    return cli


def main():
    """
    Main entry point for the CLI application.
    """
    grz_mode = os.environ.get("GRZ_MODE", "false").lower() in {"true", "yes", "1"}

    cli = build_cli(grz_mode=grz_mode)
    cli()


if __name__ == "__main__":
    main()
