"""
CLI module for handling command-line interface operations for GRZ administrators.
"""

import logging.config
from pathlib import Path

import click
import grz_common.cli as grzcli
from grz_cli.commands.submit import submit
from grz_common.cli.dump_config import dump_config
from grz_common.logging import setup_cli_logging

from . import get_versions
from .commands.archive import archive
from .commands.clean import clean
from .commands.cli_wrappers import encrypt, upload, validate
from .commands.consent import consent
from .commands.db.cli import db
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
    versions = get_versions()

    @click.group(
        cls=OrderedGroup,
        help="GRZ Control CLI for GRZ administrators.",
    )
    @click.version_option(
        version=versions["grzctl"],
        prog_name="grzctl",
        message="\n".join(f"{k} {'v' + v if v else 'unknown'}" for k, v in versions.items()),
    )
    @click.option("--log-file", metavar="FILE", type=str, help="Path to log file")
    @click.option(
        "--log-level",
        default="INFO",
        type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]),
        help="Set the log level (default: INFO)",
    )
    @grzcli.config_file
    def cli(config_file: tuple[Path], log_file: str | None = None, log_level: str = "INFO"):
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
    cli.add_command(dump_config)

    return cli


def main():
    """
    Main entry point for the CLI application.
    """
    cli = build_cli()
    cli()


if __name__ == "__main__":
    main()
