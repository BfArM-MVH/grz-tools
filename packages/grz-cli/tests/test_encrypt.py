"""Tests for ``grz-cli encrypt``, which loads the submitter private key with its configured passphrase."""

from pathlib import Path

import click.testing
import crypt4gh.keys
import crypt4gh.keys.c4gh
import pytest
import yaml
from grz_cli.cli import build_cli
from pytest_mock import MockerFixture

PASSPHRASE = "right-passphrase"


@pytest.fixture
def no_prompt(monkeypatch):
    """Fail the test if the passphrase prompt opens, and set a wrong ``C4GH_PASSPHRASE``.

    The configured passphrase comes first, so the wrong one in the environment must not matter.
    """
    monkeypatch.setenv("C4GH_PASSPHRASE", "wrong-passphrase")

    def _fail(*args, **kwargs):
        raise AssertionError("the passphrase prompt must not open")

    monkeypatch.setattr("grz_common.utils.crypt.getpass", _fail)


@pytest.mark.parametrize("inline", [False, True], ids=["path", "inline"])
def test_encrypt_signs_with_the_submitter_key_and_its_configured_passphrase(
    tmp_path: Path, submission_dir: Path, mocker: MockerFixture, no_prompt, inline: bool
):
    crypt4gh.keys.c4gh.generate(tmp_path / "grz.sec", tmp_path / "grz.pub", None, comment=None)
    private_key_path = tmp_path / "le.sec"
    crypt4gh.keys.c4gh.generate(private_key_path, tmp_path / "le.pub", PASSPHRASE.encode(), comment=None)
    keys = {"grz_public_key_path": str(tmp_path / "grz.pub"), "submitter_private_key_passphrase": PASSPHRASE}
    if inline:
        keys["submitter_private_key"] = private_key_path.read_text()
    else:
        keys["submitter_private_key_path"] = str(private_key_path)
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.dump({"keys": keys}))
    worker_encrypt = mocker.patch("grz_cli.commands.encrypt.Worker.encrypt")

    result = click.testing.CliRunner().invoke(
        build_cli(),
        ["encrypt", "--submission-dir", str(submission_dir), "--config-file", str(config_path)],
        catch_exceptions=False,
    )

    assert result.exit_code == 0, result.output
    signing_key = worker_encrypt.call_args.kwargs["submitter_private_key"]
    assert signing_key.private_bytes_raw() == crypt4gh.keys.get_private_key(private_key_path, lambda: PASSPHRASE)
