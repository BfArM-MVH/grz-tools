import logging
import sys
from datetime import UTC, datetime
from importlib.metadata import version

from grz_common.models.s3 import S3Options
from grz_common.models.version import VersionFile, VersionInfo
from packaging import version as pkg_version

logger = logging.getLogger(__name__)


def _select_active_policy(
    policies: list[VersionInfo],
    now: datetime,
) -> VersionInfo | None:
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
        return None

    applicable = [p for p in policies if p.enforced_from <= now]

    if not applicable:
        return None

    return max(applicable, key=lambda p: p.enforced_from)


def check_version_and_exit_if_needed(
    s3_options: S3Options,
    version_file_key: str = "version.json",
) -> None:
    """Validate the installed grz-cli version against the policy defined in version.json."""
    version_file = VersionFile.from_s3(s3_options, version_file_key)

    current_version = pkg_version.Version(version("grz-cli"))
    _check_policy_and_exit_if_needed("grz-cli", current_version, version_file.grzcli_version)


def check_metadata_version_and_exit_if_needed(
    s3_options: S3Options,
    metadata_schema_version: str,
    version_file_key: str = "version.json",
) -> None:
    """Validate the metadata schema version against the policy defined in version.json."""
    version_file = VersionFile.from_s3(s3_options, version_file_key)

    current_version = pkg_version.Version(metadata_schema_version)
    _check_policy_and_exit_if_needed("metadata", current_version, version_file.metadata_version)


def _check_policy_and_exit_if_needed(
    subject: str,
    current_version: pkg_version.Version,
    policies: list[VersionInfo],
) -> None:
    """Validate a version against the active policy for a subject."""
    now = datetime.now(UTC)

    policy = _select_active_policy(policies, now)
    if policy is None:
        logger.debug("No active version policy found — skipping version check.")
        return

    minimal_version = policy.minimal_version
    recommended_version = policy.recommended_version
    max_version = policy.max_version
    enforced_from = policy.enforced_from

    logger.debug(f"Current {subject} version: {current_version}")
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
            f"Your {subject} version ({current_version}) is not supported. Minimum required version is {minimal_version}."
        )
        sys.exit(1)

    # supported but behind recommended — skip if recommended_version not set
    if recommended_version is not None and minimal_version <= current_version < recommended_version:
        logger.warning(
            f"You are using {subject} {current_version}, but the recommended version is "
            f"{recommended_version}. Upgrading is strongly recommended."
        )
        return

    # too new — skip if max_version not set
    if max_version is not None and current_version > max_version:
        logger.error(
            f"{subject} version {current_version} is newer than the maximum supported version ({max_version})."
        )
        sys.exit(1)

    # best case
    logger.info(f"{subject} {current_version} is within the supported and tested range.")
