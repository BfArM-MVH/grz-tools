from .v1 import *  # noqa: F403


def get_accepted_versions() -> set[str]:
    return {"1.2.1", "1.2.2", "1.3", "1.3.0"}
