"""Shared fixtures for grz-cli tests."""

from pathlib import Path

import pytest
import yaml

MINIMAL_S3_CONFIG: dict[str, dict[str, str]] = {
    "s3": {
        "endpoint_url": "https://example.invalid",
        "bucket": "test-bucket",
        "access_key": "test-key",
        "secret": "test-secret",
    },
}


@pytest.fixture
def s3_config_path(tmp_path: Path) -> Path:
    """Write a minimal S3 config file and return its path."""
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.dump(MINIMAL_S3_CONFIG))
    return config_path


@pytest.fixture
def submission_dir(tmp_path: Path) -> Path:
    """Empty submission directory with the four standard subdirectories."""
    submission_dir = tmp_path / "submission"
    for sub in ("metadata", "files", "encrypted_files", "logs"):
        (submission_dir / sub).mkdir(parents=True)
    return submission_dir
