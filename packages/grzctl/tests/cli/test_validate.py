"""Tests for the grzctl ``validate`` wrapper command."""

import click.testing
import grz_cli.commands.validate as grz_cli_validate
import grzctl.cli
import pytest


@pytest.fixture
def grzctl_config_path(tmp_path):
    import yaml

    config = {
        "s3": {"inboxes": {"000000000": {"inbox": {"private_key_path": "/dev/null"}}}},
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
    ("flag", "expected"),
    [
        ("--mmap", True),
        ("--no-mmap", False),
        (None, False),  # default is no-mmap
    ],
)
def test_validate_forwards_mmap_to_inner_callback(tmp_path, monkeypatch, flag, expected, grzctl_config_path):
    """The grzctl ``validate`` wrapper must forward the ``mmap`` flag to the inner
    grz-cli ``validate`` callback under the parameter name it actually expects
    (``mmap``), not ``no_mmap``.

    Regression test for ``TypeError: validate() missing 1 required positional
    argument: 'mmap'``: the wrapper invokes the inner callback directly (bypassing
    Click's option parsing), so a wrong keyword name silently lands in ``**kwargs``
    and leaves the required ``mmap`` parameter unbound.
    """
    received: dict[str, object] = {}

    def fake_callback(*, configuration, submission_dir, force, threads, mmap, **kwargs):
        # Mirror the real inner callback's signature so that passing a wrong keyword
        # (e.g. ``no_mmap``) reproduces the same TypeError observed in production.
        received["mmap"] = mmap

    monkeypatch.setattr(grz_cli_validate.validate, "callback", fake_callback)

    args = ["--config", str(grzctl_config_path), "validate", "--submission-dir", str(tmp_path), "--no-update-db"]
    if flag is not None:
        args.append(flag)

    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()
    result = runner.invoke(cli, args)

    assert result.exit_code == 0, result.output
    assert received["mmap"] is expected
