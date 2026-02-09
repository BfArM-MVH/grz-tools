"""Utility functions for redacting sensitive information from files."""

from __future__ import annotations

import re
from pathlib import Path


def redact_file_patterns(file_path: Path, patterns: list[tuple[str, str]]) -> bool:
    """
    Redact patterns from a file in-place.

    :param file_path: Path to file to redact
    :param patterns: List of (pattern, replacement) tuples
    :returns: True if file was modified, False otherwise
    """
    if not patterns:
        return False

    try:
        content = file_path.read_text()
        modified = False

        for pattern_str, replacement in patterns:
            pattern = re.compile(re.escape(pattern_str))
            if pattern.search(content):
                content = pattern.sub(replacement, content)
                modified = True

        if modified:
            file_path.write_text(content)

        return modified

    except Exception:
        # let caller handle logging
        raise
