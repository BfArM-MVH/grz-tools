"""The example configs in the package docs load through the models that read them."""

from pathlib import Path
from typing import Any

import yaml
from grz_cli.models.config import EncryptConfig, UploadConfig, ValidateConfig
from grzctl.models.config import GrzctlConfig

PACKAGES_DIR = Path(__file__).parent.parent / "packages"
PLACEHOLDER_DIR = "/path/to/"


def _with_existing_files(node: Any, key_dir: Path) -> Any:
    """Replace each ``/path/to/`` placeholder with a new empty file, because path fields must name an existing file.

    :param node: Part of an example config.
    :param key_dir: Directory for the new files.
    :returns: A copy of *node* with the replaced paths.
    """
    if isinstance(node, dict):
        return {key: _with_existing_files(value, key_dir) for key, value in node.items()}
    if isinstance(node, str) and node.startswith(PLACEHOLDER_DIR):
        path = key_dir / node.removeprefix(PLACEHOLDER_DIR)
        path.touch()
        return str(path)
    return node


def _load_example(package: str, key_dir: Path) -> dict[str, Any]:
    """Read the example config of a package.

    :param package: Directory name of the package.
    :param key_dir: Directory for the files that the example's paths name.
    :returns: The example config, with paths to existing files.
    """
    example = yaml.safe_load((PACKAGES_DIR / package / "docs" / "example-config.yaml").read_text())
    return _with_existing_files(example, key_dir)


def _leaf_paths(node: dict[str, Any], prefix: tuple[str, ...] = ()) -> set[tuple[str, ...]]:
    """Collect the key path of every value that is not a mapping.

    :param node: A config as a dict.
    :param prefix: Key path of *node*.
    :returns: The key paths.
    """
    paths: set[tuple[str, ...]] = set()
    for key, value in node.items():
        if isinstance(value, dict):
            paths |= _leaf_paths(value, (*prefix, key))
        else:
            paths.add((*prefix, key))
    return paths


def test_grz_cli_example_config_loads(tmp_path: Path):
    """The models ignore unknown keys, so only the comparison catches an example key that a model no longer has."""
    example = _load_example("grz-cli", tmp_path)

    loaded: dict[str, Any] = {}
    for model in (ValidateConfig, EncryptConfig, UploadConfig):
        loaded |= model.model_validate(example).model_dump(mode="json")

    assert _leaf_paths(example) <= _leaf_paths(loaded)


def test_grzctl_example_config_loads(tmp_path: Path):
    """The models ignore unknown keys, so only the comparison catches an example key that a model no longer has."""
    example = _load_example("grzctl", tmp_path)

    loaded = GrzctlConfig.from_configuration(example).model_dump(mode="json")

    assert _leaf_paths(example) <= _leaf_paths(loaded)
