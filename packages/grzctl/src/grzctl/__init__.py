"""
GRZ Control CLI for GRZ administrators.
"""

from importlib.metadata import PackageNotFoundError, version

__version__ = "2.2.0"  # This version is managed by release-please


def get_versions() -> dict[str, str | None]:
    def get_version(package_name: str) -> str | None:
        try:
            return version(package_name)
        except PackageNotFoundError:
            return None

    return {
        "grzctl": get_version("grzctl"),
        "grz-cli": get_version("grz-cli"),
        "grz-common": get_version("grz-common"),
        "grz-db": get_version("grz-db"),
        "grz-pydantic-models": get_version("grz-pydantic-models"),
        "grz-check": get_version("grz-check"),
    }
