from .v1 import *  # noqa: F403


def get_supported_versions() -> set[str]:
    return {"1.1.1", "1.1.4"}
