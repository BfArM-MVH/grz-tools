import logging
import sys
from importlib.metadata import version

import click
from grz_common.models.s3 import S3Options
from grz_common.transfer import get_version_info
from packaging import version as pkg_version

logger = logging.getLogger(__name__)


def check_version_and_exit_if_needed(s3_options: S3Options, version_file_path: str = "version.json") -> None:
    """
    Check grz-cli version against the requirements in the version file.
    """
    # Fetch version information from S3
    version_info = get_version_info(s3_options, version_file_path)
    current_version = pkg_version.Version(version("grz-cli"))

    logger.debug(f"Current grz-cli version: {current_version}")
    logger.debug(
        f"Version file: minimal={version_info.minimal_version}, recommended={version_info.recommended_version}"
    )

    # case when the version file exists but contains outdated version info and GRZ needs to fix this.
    if version_info.recommended_version < current_version:
        msg = (
            f"The version file in S3 appears outdated — it lists recommended_version={version_info.recommended_version}, "
            f"but you are running grz-cli {current_version}. "
            "This means GRZ needs to update the version file to reflect the recommended CLI versions."
        )
        logger.critical(msg)
        sys.exit(1)

    # case when the CLI is too old and the LE must upgrade and possibly revalidate metadata.
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
        sys.exit(1)

    # case when the version is behind the latest but still supported
    elif version_info.minimal_version <= current_version < version_info.recommended_version:
        msg = (
            f"You are using grz-cli {current_version}, while the recommended version is {version_info.recommended_version}.\n"
            "It is recommended to upgrade to the latest version for the newest features and bug fixes."
        )
        logger.warning(msg)

    # version is up-to-date
    else:
        msg = f"grz-cli {current_version} is up to date."
        logger.info(msg)
