import logging
import sys
from datetime import UTC, date, datetime
from importlib.metadata import version

from grz_common.models.s3 import S3Options
from grz_common.models.version import VersionFile, VersionInfo
from packaging import version as pkg_version

logger = logging.getLogger(__name__)

def _select_active_policy(
    policies: list[VersionInfo],
    today: date,
) -> VersionInfo:
    """Select the version policy that is active for the given date.

    :param policies: A non-empty list of available version policies.
                     Each policy defines version constraints and an
                     enforcement start date.
    :param today: The reference date used to determine which policy
                  is currently active (typically ``date.today()``).
    :returns: The policy that should be applied for the given date.
    :raises ValueError: If policies is empty.
    """
    if not policies:
        raise ValueError("No version policies defined.")

    applicable = [p for p in policies if p.enforced_from <= today]

    if not applicable:
        # No policy has become active yet, so return the earliest one
        return min(policies, key=lambda p: p.enforced_from)

    # Return the most recent policy that is already active
    return max(applicable, key=lambda p: p.enforced_from)


def check_version_and_exit_if_needed(
    s3_options: S3Options,
    version_file_key: str = "version.json",
) -> None:
    """Validate the installed grz-cli version against the policy defined in version.json."""
    version_file = VersionFile.from_s3(s3_options, version_file_key)

    current_version = pkg_version.Version(version("grz-cli"))
    today = datetime.now(UTC).date()

    policy = _select_active_policy(version_file.grzcli_version, today)

    minimal_version = policy.minimal_version
    recommended_version = policy.recommended_version
    max_version = policy.max_version
    enforced_from = policy.enforced_from

    logger.debug(f"Current grz-cli version: {current_version}")
    logger.debug(
        f"Active policy: "
        f"minimal={minimal_version}, "
        f"recommended={recommended_version}, "
        f"max_version={max_version}, "
        f"enforced_from={enforced_from}"
    )

    # old version
    if current_version < minimal_version:
        logger.error(
            f"Your grz-cli version ({current_version}) is not supported. "
            f"Minimum required version is {minimal_version}."
        )
        sys.exit(1)

    # supported but behind recommended
    if minimal_version <= current_version < recommended_version:
        logger.warning(
            f"You are using grz-cli {current_version}, but the recommended version is "
            f"{recommended_version}. Upgrading is strongly recommended."
        )
        return

    # too new
    if current_version > max_version:
        logger.error(
            f"grz-cli version {current_version} is newer than the maximum supported "
            f"version ({max_version})."
        )
        sys.exit(1)

    # best case
    logger.info(
        f"grz-cli {current_version} is within the supported and tested range."
    )
