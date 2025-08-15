"""
CLI module for handling command-line interface operations for GRZ administrators.
"""

import logging
import logging.config
import shutil
import subprocess
from importlib.metadata import version
from textwrap import dedent

import click
from grz_cli.commands.encrypt import encrypt
from grz_cli.commands.submit import submit
from grz_cli.commands.upload import upload
from grz_cli.commands.validate import validate
from grz_common.logging import setup_cli_logging

from .commands.archive import archive
from .commands.clean import clean
from .commands.consent import consent
from .commands.db import db
from .commands.decrypt import decrypt
from .commands.download import download
from .commands.list_submissions import list_submissions
from .commands.pruefbericht import pruefbericht
from .commands.report import report

log = logging.getLogger(__name__)


class OrderedGroup(click.Group):
    """
    A click Group that keeps track of the order in which commands are added.
    """

    def list_commands(self, ctx):
        """Return the list of commands in the order they were added."""
        return list(self.commands.keys())


def build_cli():
    """
    Factory for building the CLI application.
    """

    @click.group(
        cls=OrderedGroup,
        help="GRZ Control CLI for GRZ administrators.",
    )
    @click.version_option(
        version=version("grzctl"),
        prog_name="grzctl",
        message=dedent(f"""\
        %(prog)s v%(version)s
        grz-cli v{version("grz-cli")}
        grz-common v{version("grz-common")}
        grz-db v{version("grz-db")}
        grz-pydantic-models v{version("grz-pydantic-models")}
        """)
        + (
            subprocess.run(["grz-check", "--version"], capture_output=True, text=True).stdout.strip()  # noqa: S603, S607
            if shutil.which("grz-check") is not None
            else ""
        ),
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
        setup_cli_logging(log_file, log_level)

    # For convenience, include grz-cli commands as well.
    cli.add_command(validate)
    cli.add_command(encrypt)
    cli.add_command(upload)
    cli.add_command(submit)

    cli.add_command(list_submissions, name="list")
    cli.add_command(download)
    cli.add_command(decrypt)
    cli.add_command(archive)
    cli.add_command(clean)
    cli.add_command(consent)
    cli.add_command(pruefbericht)
    cli.add_command(db)
    cli.add_command(report)

    return cli


def main():
    """
    Main entry point for the CLI application.
    """
    cli = build_cli()
    cli()


if __name__ == "__main__":
    main()
