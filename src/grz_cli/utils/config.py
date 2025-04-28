"""Configuration utilities for the GRZ CLI."""

import logging
from os import PathLike

import yaml

from ..models.config import ConfigModel

log = logging.getLogger(__name__)


def read_config(config_path: str | PathLike) -> ConfigModel:
    """Reads the configuration file and validates it against the schema."""
    with open(config_path, encoding="utf-8") as f:
        config = ConfigModel(**yaml.safe_load(f))

    return config
