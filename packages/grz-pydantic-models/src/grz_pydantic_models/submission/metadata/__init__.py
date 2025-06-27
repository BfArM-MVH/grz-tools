from .v1 import *  # noqa: F403


def get_accepted_versions() -> set[str]:
    return {"1.1.7", "1.1.8", "1.1.9"}
