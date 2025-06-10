from .std import deprecated
from .submission.metadata import get_supported_versions as _get_supported_versions

__version__ = "1.4.0"


@deprecated(
    msg="Using get_supported_versions() directly from the top level package is deprecated. Please import it from submission.metadata instead."
)
def get_supported_versions() -> set[str]:
    return _get_supported_versions()
