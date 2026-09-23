"""Tests for the grzctl ``encrypt`` command."""

from pathlib import Path
from unittest.mock import patch

import click.testing
import grzctl.cli
import pytest
import yaml


@pytest.fixture
def grzctl_config_path(tmp_path, crypt4gh_public_key):
    config = {
        "leistungserbringer": {"000000000": {"inbox_buckets": {"inbox": {"private_key_path": "/dev/null"}}}},
        "archives": {
            "consented": {"s3": {"bucket": "consented"}, "public_key": crypt4gh_public_key},
            "non_consented": {"s3": {"bucket": "non_consented"}, "public_key_path": "/dev/null"},
        },
        "db": {"database_url": "sqlite:///:memory:", "author": {"name": "test"}},
        "pruefbericht": {},
        "keys": {"grz_private_key_path": "/dev/null"},
        "identifiers": {"grz": "GRZT00000"},
    }
    config_path = tmp_path / "config.yaml"
    with open(config_path, "w") as f:
        yaml.dump(config, f)
    return config_path


def test_encrypt_uses_an_inline_archive_public_key(tmp_path, grzctl_config_path, crypt4gh_public_key):
    """The consented archive's public key may be given inline instead of as a file path.

    ``grzctl encrypt`` writes it to a temporary file, since ``Worker.encrypt`` only takes a path.
    The file is cleaned up again once ``Worker.encrypt`` returns, so its content has to be read
    from inside the mocked call rather than after ``invoke`` comes back.
    """
    submission_dir = tmp_path / "submission"
    for sub in ("metadata", "files", "logs", "encrypted_files"):
        (submission_dir / sub).mkdir(parents=True)

    used_key_content = None

    def _capture_key_content(*, recipient_public_key_path, **kwargs):
        nonlocal used_key_content
        used_key_content = Path(recipient_public_key_path).read_text()

    with patch("grzctl.commands.encrypt.Worker") as mock_worker_cls:
        mock_worker = mock_worker_cls.return_value
        mock_submission = mock_worker.parse_submission.return_value
        mock_submission.metadata.content.submission_id = "S1"
        mock_submission.metadata.content.consents_to_research.return_value = True
        mock_worker.encrypt.side_effect = _capture_key_content

        runner = click.testing.CliRunner()
        cli = grzctl.cli.build_cli()
        result = runner.invoke(
            cli,
            [
                "--config",
                str(grzctl_config_path),
                "encrypt",
                "--submission-dir",
                str(submission_dir),
                "--no-update-db",
            ],
        )

        assert result.exit_code == 0, result.output
        mock_worker.encrypt.assert_called_once()
        assert used_key_content == crypt4gh_public_key
