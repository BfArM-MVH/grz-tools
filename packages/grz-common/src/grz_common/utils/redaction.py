"""Utility functions for redacting sensitive information from files."""

import re
from pathlib import Path


def redact_file(input_file: Path, output_file: Path, patterns: list[tuple[str, str]]) -> bool:
    """
    Redact patterns from a file in-place.

    :param input_file: Path to file to redact
    :param output_file: Path to file to write the redacted file to
    :param patterns: List of (pattern, replacement) tuples
    :returns: True if file was modified, False otherwise
    """
    if not patterns:
        return False

    content = input_file.read_text()
    modified = False

    for pattern_str, replacement in patterns:
        pattern = re.compile(re.escape(pattern_str))
        if pattern.search(content):
            content = pattern.sub(replacement, content)
            modified = True

    output_file.write_text(content)

    return modified
