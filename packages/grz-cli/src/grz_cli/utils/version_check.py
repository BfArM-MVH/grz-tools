import logging
import sys
from datetime import UTC, datetime  # removed date
from importlib.metadata import version

from grz_common.models.s3 import S3Options
from grz_common.models.version import VersionFile, VersionInfo
from packaging import version as pkg_version

logger = logging.getLogger(__name__)


def _select_active_policy(
    policies: list[VersionInfo],
    now: datetime,  # changed
) -> VersionInfo:
    """Select the version policy that is active for the given datetime.

    :param policies: A non-empty list of available version policies.
                     Each policy defines version constraints and an
                     enforcement start datetime.
    :param now: The reference datetime used to determine which policy
                is currently active (typically ``datetime.now(UTC)``).
    :returns: The policy that should be applied for the given datetime.
    :raises ValueError: If policies is empty.
    """
    if not policies:
        raise ValueError("No version policies defined.")

    applicable = [p for p in policies if p.enforced_from <= now]  # changed

    if not applicable:
        return min(policies, key=lambda p: p.enforced_from)

    return max(applicable, key=lambda p: p.enforced_from)


def check_version_and_exit_if_needed(
    s3_options: S3Options,
    version_file_key: str = "version.json",
) -> None:
    """Validate the installed grz-cli version against the policy defined in version.json."""
    version_file = VersionFile.from_s3(s3_options, version_file_key)

    current_version = pkg_version.Version(version("grz-cli"))
    now = datetime.now(UTC)  # changed

    policy = _select_active_policy(version_file.grzcli_version, now)  # changed

    if policy is None:
        logger.debug("No active version policy found — skipping version check.")
        return

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
            f"Your grz-cli version ({current_version}) is not supported. Minimum required version is {minimal_version}."
        )
        sys.exit(1)

    # supported but behind recommended — skip if recommended_version not set
    if recommended_version is not None and minimal_version <= current_version < recommended_version:  # changed
        logger.warning(
            f"You are using grz-cli {current_version}, but the recommended version is "
            f"{recommended_version}. Upgrading is strongly recommended."
        )
        return

    # too new — skip if max_version not set
    if max_version is not None and current_version > max_version:  # changed
        logger.error(f"grz-cli version {current_version} is newer than the maximum supported version ({max_version}).")
        sys.exit(1)

    # best case
    logger.info(f"grz-cli {current_version} is within the supported and tested range.")
