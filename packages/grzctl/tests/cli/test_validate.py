"""Tests for the grzctl ``validate`` command."""

from unittest.mock import patch

import click.testing
import grzctl.cli
import pytest


@pytest.fixture
def grzctl_config_path(tmp_path):
    import yaml

    config = {
        "leistungserbringer": {"000000000": {"inbox_buckets": {"inbox": {"private_key_path": "/dev/null"}}}},
        "archives": {
            "consented": {"s3": {"bucket": "consented"}, "public_key_path": "/dev/null"},
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


@pytest.mark.parametrize(
    ("flag", "expected_no_mmap"),
    [
        ("--mmap", False),
        ("--no-mmap", True),
        (None, True),  # default is no-mmap
    ],
)
def test_validate_forwards_mmap_to_worker(tmp_path, flag, expected_no_mmap, grzctl_config_path):
    """The grzctl ``validate`` command must forward the ``--mmap/--no-mmap`` flag to
    ``Worker.validate`` (as the inverted ``no_mmap`` keyword).

    Regression test for a previous wrapper implementation which invoked grz-cli's
    inner callback directly and silently dropped misnamed keywords into ``**kwargs``.
    """
    submission_dir = tmp_path / "submission"
    for sub in ("metadata", "files", "logs"):
        (submission_dir / sub).mkdir(parents=True)

    with patch("grzctl.commands.validate.Worker") as mock_worker_cls:
        mock_worker = mock_worker_cls.return_value
        mock_worker.parse_submission.return_value.metadata.content.submission_id = "S1"

        args = [
            "--config",
            str(grzctl_config_path),
            "validate",
            "--submission-dir",
            str(submission_dir),
            "--submitter-id",
            "000000000",
            "--no-update-db",
        ]
        if flag is not None:
            args.append(flag)

        runner = click.testing.CliRunner()
        cli = grzctl.cli.build_cli()
        result = runner.invoke(cli, args)

        assert result.exit_code == 0, result.output
        mock_worker.validate.assert_called_once()
        assert mock_worker.validate.call_args.kwargs["no_mmap"] is expected_no_mmap
