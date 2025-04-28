"""
CLI module for handling command-line interface operations.
"""

import importlib
import logging
import logging.config
import os
from os import sched_getaffinity
from pathlib import Path

import click
import grz_pydantic_models.submission.metadata
import platformdirs

from .commands.decrypt import decrypt
from .commands.download import download
from .commands.encrypt import encrypt
from .commands.list_submissions import list_submissions
from .commands.submit import submit
from .commands.upload import upload
from .commands.validate import validate
from .constants import PACKAGE_ROOT
from .logging_setup import add_filelogger

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


class OrderedGroup(click.Group):
    """
    A click Group that keeps track of the order in which commands are added.
    """

    def list_commands(self, ctx):
        """Return the list of commands in the order they were added."""
        return list(self.commands.keys())


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

    # Add commands
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
