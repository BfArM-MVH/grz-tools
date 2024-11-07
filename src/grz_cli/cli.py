"""
CLI module for handling command-line interface operations.
"""

import logging
import logging.config
import sys
from os import PathLike
from pathlib import Path

import click
import platformdirs
import yaml

from .constants import PACKAGE_ROOT
from .logging_setup import add_filelogger
from .models.config import ConfigModel
from .worker import Worker

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
FILE_R_E = click.Path(
    exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True
)

submission_dir = click.option(
    "-s",
    "--submission-dir",
    metavar="PATH",
    type=DIR_R_E,
    required=True,
    help="Path to the submission directory containing both 'metadata/' and 'files/' directories",
)


metadata_dir = click.option(
    "-m",
    "--metadata-dir",
    metavar="PATH",
    type=DIR_RW_C,
    required=False,
    default=None,
    help="Path to the metadata directory containing the metadata.json file",
)

files_dir = click.option(
    "-f",
    "--files-dir",
    metavar="PATH",
    type=DIR_RW_C,
    required=False,
    default=None,
    help="Path to the files linked in the submission",
)

working_dir = click.option(
    "-w",
    "--working-dir",
    metavar="PATH",
    type=DIR_RW_E,
    required=False,
    default=None,
    callback=lambda c, p, v: v or c.params.get("submission_dir"),
    help="Path to a working directory where intermediate files can be stored",
)

config_file = click.option(
    "-c",
    "--config-file",
    metavar="STRING",
    type=FILE_R_E,
    required=True,
    default=DEFAULT_CONFIG_PATH,
    help="Path to config file",
)

encrypted_files_dir = click.option(
    "-o",
    "--encrypted-files-dir",
    metavar="PATH",
    type=DIR_RW_C,
    required=False,
    default=None,
    help="Path to a directory where the encrypted files can be stored",
)

decrypted_files_dir = click.option(
    "-o",
    "--decrypted-files-dir",
    metavar="STRING",
    type=DIR_RW_C,
    required=False,
    default=None,
    help="Path to a directory where the decrypted files can be stored",
)

threads = click.option(
    "--threads",
    default=None,
    type=int,
    help="Number of threads to use for parallel operations",
)


class OrderedGroup(click.Group):
    """
    A click Group that keeps track of the order in which commands are added.
    """

    def list_commands(self, ctx):
        """Return the list of commands in the order they were added."""
        return list(self.commands.keys())


@click.group(
    cls=OrderedGroup,
    help="Validate, encrypt, decrypt and upload submissions to a GRZ/GDC.",
)
@click.version_option(version="0.1", prog_name="grz-cli")
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


@cli.command()
@submission_dir
@metadata_dir
@files_dir
@working_dir
@threads
def validate(
    submission_dir: str,
    metadata_dir: str,
    files_dir: str,
    working_dir: str,
    threads: int | None,
):
    """
    Validates the sha256 checksum of the sequence data files. This command must be executed
    before the encryption and upload can start.
    """
    log.info("Starting validation...")

    working_dir_path: Path = Path(working_dir)
    submission_dir_path: Path = Path(submission_dir)

    worker_inst = Worker(
        working_dir=working_dir_path,
        files_dir=submission_dir_path / "files" if files_dir is None else files_dir,
        metadata_dir=(metadata_dir or submission_dir_path / "metadata"),
        threads=threads,
    )
    worker_inst.validate()

    log.info("Validation finished!")


@cli.command()
@submission_dir
@metadata_dir
@files_dir
@working_dir
@encrypted_files_dir
@config_file
@threads
def encrypt(  # noqa: PLR0913
    submission_dir,
    working_dir,
    metadata_dir,
    files_dir,
    encrypted_files_dir,
    config_file,
    threads,
):
    """
    Encrypt a submission using the GRZ public key.
    """
    config = read_config(config_file)

    grz_pubkey_path = config.grz_public_key_path

    submitter_privkey_path = config.submitter_private_key_path
    if submitter_privkey_path == "":
        submitter_privkey_path = None

    log.info("Starting encryption...")

    working_dir = Path(working_dir)
    submission_dir = Path(submission_dir)

    worker_inst = Worker(
        working_dir=working_dir,
        encrypted_files_dir=encrypted_files_dir,
        files_dir=submission_dir / "files" if files_dir is None else files_dir,
        metadata_dir=(
            submission_dir / "metadata" if metadata_dir is None else metadata_dir
        ),
        threads=threads,
    )
    worker_inst.encrypt(
        grz_pubkey_path, submitter_private_key_path=submitter_privkey_path
    )

    log.info("Encryption successful!")


@cli.command()
@submission_dir
@metadata_dir
@encrypted_files_dir
@working_dir
@decrypted_files_dir
@config_file
@threads
def decrypt(  # noqa: PLR0913
    submission_dir: str,
    metadata_dir: str,
    encrypted_files_dir: str,
    working_dir: str,
    decrypted_files_dir: str,
    config_file: str,
    threads: int | None,
):
    """
    Decrypt a submission using the GRZ private key.
    """
    config = read_config(config_file)

    grz_privkey_path = config.grz_private_key_path
    if not grz_privkey_path:
        log.error("GRZ private key path is required for decryption.")
        sys.exit(1)

    log.info("Starting encryption...")

    working_dir_path: Path = Path(working_dir)
    submission_dir_path: Path = Path(submission_dir)

    worker_inst = Worker(
        working_dir=working_dir_path,
        files_dir=(
            working_dir_path / "files"
            if decrypted_files_dir is None
            else decrypted_files_dir
        ),
        metadata_dir=(metadata_dir or submission_dir_path / "metadata"),
        encrypted_files_dir=(
            submission_dir_path / "encrypted_files"
            if encrypted_files_dir is None
            else encrypted_files_dir
        ),
        threads=threads,
    )
    worker_inst.decrypt(grz_privkey_path)

    log.info("Encryption successful!")


@cli.command()
@submission_dir
@metadata_dir
@encrypted_files_dir
@working_dir
@config_file
@threads
def upload(  # noqa: PLR0913
    submission_dir,
    metadata_dir,
    encrypted_files_dir,
    working_dir,
    config_file,
    threads,
):
    """
    Upload a submission file to s3 using the provided configuration.
    :param submission_dir: The path to the encrypted submission directory.
    :param metadata_dir: The path to the metadata directory.
    :param encrypted_files_dir: The path to the encrypted files directory.
    :param working_dir: Path to a working directory where intermediate files can be stored
    :param config_file: The path to the configuration file.
    """
    config = read_config(config_file)

    log.info("Starting upload...")
    working_dir_path = Path(working_dir)
    submission_dir_path = Path(submission_dir)

    worker_inst = Worker(
        working_dir=working_dir_path,
        metadata_dir=(metadata_dir or submission_dir_path / "metadata"),
        encrypted_files_dir=(
            submission_dir_path / "encrypted_files"
            if encrypted_files_dir is None
            else encrypted_files_dir
        ),
        threads=threads,
    )
    worker_inst.upload(config)

    log.info("Upload finished!")


@cli.command()
@click.argument("submission_id", type=str)
@submission_dir
@metadata_dir
@encrypted_files_dir
@config_file
def download(
    submission_id,
    submission_dir,
    metadata_dir,
    encrypted_files_dir,
    config_file,
):
    """
    Download a submission file from s3 using the provided configuration.
    :param submission_id: Submission id
    :param submission_dir: The path to the encrypted submission directory.
    :param metadata_dir: The path to the metadata directory.
    :param encrypted_files_dir: The path to the encrypted files directory.
    :param config_file: The path to the configuration file.
    """
    config = read_config(config_file)

    log.info("Starting download...")

    submission_dir_path = Path(submission_dir)

    worker_inst = Worker(
        metadata_dir=(
            submission_dir_path / submission_id / "metadata"
            if metadata_dir is None
            else metadata_dir
        ),
        encrypted_files_dir=(
            submission_dir_path / submission_id / "files"
            if encrypted_files_dir is None
            else encrypted_files_dir
        ),
    )
    worker_inst.download(config)

    log.info("Download finished!")


def read_config(config_path: str | PathLike) -> ConfigModel:
    """Reads the configuration file and validates it against the schema."""
    with open(config_path, encoding="utf-8") as f:
        config = ConfigModel(**yaml.safe_load(f))

    return config


if __name__ == "__main__":
    cli()
