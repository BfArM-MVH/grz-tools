import logging
import sys
from importlib.metadata import version

from grz_common.models.s3 import S3Options
from grz_common.models.version import VersionFile
from packaging import version as pkg_version

logger = logging.getLogger(__name__)


def check_version_and_exit_if_needed(s3_options: S3Options, version_file_path: str = "version.json") -> None:
    """
    Check grz-cli version against the requirements in the version file.
    """
    # Fetch version information from S3
    version_info = VersionFile.from_s3(s3_options, version_file_path)

    current_version = pkg_version.Version(version("grz-cli"))
    minimal_version = pkg_version.Version(version_info.minimal_version)
    recommended_version = pkg_version.Version(version_info.recommended_version)

    logger.debug(f"Current grz-cli version: {current_version}")
    logger.debug(
        f"Version file: minimal={minimal_version}, recommended={recommended_version}"
    )

    # case when the version file exists but contains outdated version info
    if recommended_version < current_version:
        msg = (
            f"The version file in S3 appears outdated — it lists recommended_version={recommended_version}, "
            f"but you are running grz-cli {current_version}. "
            "This means GRZ needs to update the version file to reflect the recommended CLI versions."
        )
        logger.critical(msg)
        sys.exit(1)

    # case when the CLI is too old
    elif current_version < minimal_version:
        msg = (
            f"Your grz-cli version ({current_version}) is too old and not supported anymore.\n"
            f"Minimum required version is {minimal_version}.\n\n"
            "Please upgrade grz-cli to the latest version before proceeding.\n"
            "After upgrading, you may need to:\n"
            "  - Update your metadata files to the new schema\n"
            "  - Re-run validation before uploading again\n"
            "  - Check release notes for any breaking changes\n"
        )
        logger.error(msg)
        sys.exit(1)

    # version is behind recommended but still supported
    elif minimal_version <= current_version < recommended_version:
        msg = (
            f"You are using grz-cli {current_version}, while the recommended version is {recommended_version}.\n"
            "It is recommended to upgrade to the latest version for the newest features and bug fixes."
        )
        logger.warning(msg)

    # version is up-to-date
    else:
        msg = f"grz-cli {current_version} is up to date."
        logger.info(msg)

