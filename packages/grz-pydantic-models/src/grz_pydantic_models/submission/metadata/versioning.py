"""
versioning

Handles aspects related to submission metadata schema versioning, such as comparison.
"""

from functools import total_ordering
from itertools import zip_longest
from typing import Any


@total_ordering
class Version:
    """
    Represents a version of the GRZ submission metadata.

    Currently, they can only be simple versions like 1, 1.2, 1.2.1, etc.
    """

    def __init__(self, version: str):
        try:
            self._components = map(int, version.split("."))
        except Exception as err:
            raise ValueError(f"Failed to parse '{version}' as a version string") from err

    def __eq__(self, other: Any) -> bool:
        """Equal to"""
        if not isinstance(other, Version):
            return False
        return all((a == b) for a, b in zip_longest(self._components, other._components, fillvalue=0))

    def __lt__(self, other: Any) -> bool:
        """Less than"""
        if not isinstance(other, Version):
            raise TypeError(f"Version can only be compared to another Version, not {type(other)}")
        for a, b in zip_longest(self._components, other._components, fillvalue=0):
            if a != b:
                return a < b
        return False
