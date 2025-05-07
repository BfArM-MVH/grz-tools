"""
CLI module for handling command-line interface operations.
"""

import logging
import logging.config
import os
from importlib.metadata import version

import click
import grz_pydantic_models.submission.metadata

from .commands.archive import archive
from .commands.clean import clean
from .commands.consent import consent
from .commands.decrypt import decrypt
from .commands.download import download
from .commands.encrypt import encrypt
from .commands.list_submissions import list_submissions
from .commands.pruefbericht import pruefbericht
from .commands.submit import submit
from .commands.upload import upload
from .commands.validate import validate
from .constants import PACKAGE_ROOT
from .logging_setup import add_filelogger

log = logging.getLogger(PACKAGE_ROOT + ".cli")


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
        version=version("grz-cli"),
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
        cli.add_command(archive)
        cli.add_command(clean)
        cli.add_command(consent)
        cli.add_command(pruefbericht)

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
