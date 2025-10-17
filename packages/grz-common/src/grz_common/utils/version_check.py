import sys
from importlib.metadata import version
import logging
import click

from grz_common.models.s3 import S3Options
from grz_common.transfer import get_version_info

logger = logging.getLogger(__name__)

def check_version_and_exit_if_needed(s3_options: S3Options, version_file_path: str = "version.json") -> None:
    """
    Check grz-cli version against the requirements in the version file.

    Behavior:
    case 1 - the version file is missing or inaccessible, `get_version_info` raises a hard error.
    case 2 - if the version file exists but contains outdated version info, GRZ staff must fix it.
    case 3 - if the CLI is too old, the LE must upgrade and possibly revalidate metadata.
    """

    # Fetch version information from S3
    version_info = get_version_info(s3_options, version_file_path)
    current_version = version("grz-cli")

    logger.debug(f"Current grz-cli version: {current_version}")
    logger.debug(f"Version file: minimal={version_info.minimal_version}, latest={version_info.latest_version}")

    # case 2
    if version_info.latest_version < current_version:
        msg = (
            f"The version file in S3 appears outdated — it lists latest_version={version_info.latest_version}, "
            f"but you are running grz-cli {current_version}. "
            "This means GRZ staff must update the version file to reflect the latest supported CLI versions."
        )
        logger.critical(msg)
        click.echo(click.style(f"ERROR: {msg}", fg="red"))
        sys.exit(1)

    # case 3
    elif current_version < version_info.minimal_version:
        msg = (
            f"Your grz-cli version ({current_version}) is too old and not supported anymore.\n"
            f"Minimum required version is {version_info.minimal_version}.\n\n"
            "Please upgrade grz-cli to the latest version before proceeding.\n"
            "After upgrading, you may need to:\n"
            "  - Update your metadata files to the new schema\n"
            "  - Re-run validation before uploading again\n"
            "  - Check release notes for any breaking changes\n"
        )
        logger.error(msg)
        click.echo(click.style(f"ERROR: {msg}", fg="red"))
        sys.exit(1)

    # version is behind latest but still supported
    elif version_info.minimal_version <= current_version < version_info.latest_version:
        msg = (
            f"You are using grz-cli {current_version}, while the latest version is {version_info.latest_version}.\n"
            "It is recommended to upgrade to the latest version for the newest features and bug fixes."
        )
        logger.warning(msg)
        click.echo(click.style(f"WARNING: {msg}", fg="yellow"))

    # version is up-to-date
    else:
        msg = f"grz-cli {current_version} is up to date."
        logger.info(msg)
        click.echo(click.style(msg, fg="green"))
