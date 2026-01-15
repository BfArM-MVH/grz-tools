import json
from copy import deepcopy
from pathlib import Path

import yaml

__all__ = [
    "merge_config_dicts",
    "read_and_merge_config_files",
]


def _merge_config_dicts_recursive(a: dict, b: dict, path) -> dict:
    """Helper function to merge two configuration dictionaries recursively."""
    for key, b_val in b.items():
        # If key exists in the target dictionary, merge values
        if key in a:
            a_val = a[key]

            if b_val is None:
                continue

            if a_val is None:
                a[key] = b_val
                continue

            # If both values are dictionaries, merge them recursively
            if isinstance(a_val, dict) and isinstance(b_val, dict):
                _merge_config_dicts_recursive(a_val, b_val, [*path, str(key)])
            # If the value type from ``a`` matches the value type of ``b``, ``b`` replaces the value in ``a``.
            elif type(a_val) == type(b_val):
                # Use deepcopy to avoid modifying the original objects
                a[key] = deepcopy(b_val)
            # Conflicting values
            else:
                raise ValueError(
                    "Conflict at " + ".".join([*path, str(key)]) + ": " + repr(a_val) + " != " + repr(b_val)
                )
        else:
            # If key does not exist in the target dictionary, add it
            a[key] = b_val
    return a


def merge_config_dicts(a: dict, b: dict) -> dict:
    """Merge two configuration dictionaries recursively.

    This function merges dictionary ``b`` into dictionary ``a``.

    - For keys present in both ``a`` and ``b``:
        - If both values are dictionaries, they are merged recursively.
        - If the value type from ``a`` matches the value type of ``b``, ``b`` replaces the value in ``a``.
        - If one of both values is ``None``, the non-``None`` value is used.
        - If the value types are different and cannot be merged, a ``ValueError`` is raised.
    - For keys present only in ``b``, they are added to ``a``.

    :param a: The target dictionary to merge into.
    :param b: The source dictionary to merge from.
    :return: The merged dictionary ``a``.
    :raises ValueError: If there is a conflict between values that cannot be merged.
    """
    # Create a deep copy of `a` to avoid modifying the original dictionary
    a = deepcopy(a)
    return _merge_config_dicts_recursive(a, b, path=[])


def read_and_merge_config_files(config_files: list[Path]) -> dict:
    """
    Read and merge multiple configuration files in YAML or JSON format.

    :param config_files:
    :return: Merged configuration dictionary.
    :raises RuntimeError: If there is an error reading any of the configuration files.
    """
    configuration: dict[str, object] = {}
    for curr_config_file in config_files:
        try:
            if curr_config_file.suffix.lower() == ".json":
                with open(curr_config_file) as fd:
                    curr_config = json.load(fd)
            else:
                with open(curr_config_file) as fd:
                    curr_config = yaml.safe_load(fd)

            # merge configurations
            configuration = merge_config_dicts(configuration, curr_config)
        except Exception as e:
            raise RuntimeError(f"Error reading configuration file: '{curr_config_file}'") from e

    return configuration
