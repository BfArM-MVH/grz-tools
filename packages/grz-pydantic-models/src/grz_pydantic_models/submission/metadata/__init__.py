from ...std import deprecated
from .v1 import *  # noqa: F403


@deprecated(msg="get_supported_versions() is deprecated. Use get_accepted_versions() instead.")
def get_supported_versions() -> set[str]:
    return {"1.1.1", "1.1.4"}


def get_accepted_versions() -> set[str]:
    return {"1.1.7", "1.1.8"}
