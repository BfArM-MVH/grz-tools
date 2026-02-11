import logging
import sys
from datetime import UTC, datetime
from importlib.metadata import version

from grz_common.models.s3 import S3Options
from grz_common.models.version import VersionFile
from packaging import version as pkg_version

logger = logging.getLogger(__name__)


def check_version_and_exit_if_needed(
    s3_options: S3Options,
    version_file_key: str = "version.json",
) -> None:
    """Validate the installed grz-cli version against the policy defined in version.json.
    Policy components:
      - Too old (before enforcement date): warning
      - Too old (after enforcement date): exit
      - Supported but behind recommended: warning
      - Within tested range: info
      - More recent than tested: warning
    """
    version_info = VersionFile.from_s3(s3_options, version_file_key)

    current_version = pkg_version.Version(version("grz-cli"))
    schema_version = version_info.schema_version
    minimal_version = version_info.minimal_version
    recommended_version = version_info.recommended_version
    max_version = version_info.max_version
    enforced_from = version_info.enforced_from

    today = datetime.now(UTC).date()

    logger.debug(f"Current grz-cli version: {current_version}")
    logger.debug(
        f"Version policy: {schema_version}"
        f"minimal={minimal_version}, "
        f"recommended={recommended_version}, "
        f"max_tested={max_version}, "
        f"enforced_from={enforced_from}"
    )

    # really old version
    if current_version < minimal_version:
        if today < enforced_from:
            logger.warning(
                f"Your grz-cli version ({current_version}) will become unsupported on {enforced_from}. "
                f"The minimum required version will be {minimal_version}. Please upgrade soon."
            )
            return
        else:
            logger.error(
                f"Your grz-cli version ({current_version}) is not supported. "
                f"Minimum required version is {minimal_version}."
            )
            sys.exit(1)

    # supported but not the most recent
    if minimal_version <= current_version < recommended_version:
        logger.warning(
            f"You are using grz-cli {current_version}, but the recommended version is "
            f"{recommended_version}. Upgrading is strongly recommended for latest fixes and features."
        )
        return

    # too recent that has not been tested fully
    if current_version > max_version:
        logger.warning(
            f"You are running grz-cli {current_version}, which is newer than the latest "
            f"tested version ({max_version}). Note that this version has not been tested thoroughly."
        )
        return

    # best case scenario
    logger.info(f"grz-cli {current_version} is within the supported and tested range.")
