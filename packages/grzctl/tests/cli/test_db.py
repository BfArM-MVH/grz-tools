"""
Tests for grzctl db subcommand
"""

import importlib.resources
import sqlite3
from pathlib import Path

import click.testing
import cryptography.hazmat.primitives.serialization as cryptser
import grzctl.cli
import pytest
import yaml
from cryptography.hazmat.primitives.asymmetric.ed25519 import Ed25519PrivateKey
from grz_common.constants import REDACTED_TAN
from grz_pydantic_models.submission.metadata import GrzSubmissionMetadata
from grzctl.models.config import DbConfig

from .. import resources as test_resources


@pytest.fixture
def blank_database_config(tmp_path: Path) -> DbConfig:
    private_key = Ed25519PrivateKey.generate()
    private_key_path = tmp_path / "alice.sec"
    with open(private_key_path, "wb") as private_key_file:
        private_key_file.write(
            private_key.private_bytes(
                encoding=cryptser.Encoding.PEM,
                format=cryptser.PrivateFormat.OpenSSH,
                encryption_algorithm=cryptser.NoEncryption(),
            )
        )

    public_key = private_key.public_key()
    public_key_path = tmp_path / "alice.pub"
    with open(public_key_path, "wb") as public_key_file:
        public_key_file.write(
            public_key.public_bytes(encoding=cryptser.Encoding.OpenSSH, format=cryptser.PublicFormat.OpenSSH)
        )
        # add the comment too
        public_key_file.write(b" alice")

    return DbConfig(
        db={
            "database_url": "sqlite:///" + str((tmp_path / "submission.db.sqlite").resolve()),
            "author": {
                "name": "alice",
                "private_key_path": str(private_key_path.resolve()),
                "private_key_passphrase": "",
            },
            "known_public_keys": str(public_key_path.resolve()),
        }
    )


@pytest.fixture
def blank_initial_database_config_path(tmp_path: Path, blank_database_config: DbConfig) -> Path:
    config_path = tmp_path / "config.db.yaml"
    with open(config_path, "w") as config_file:
        config_file.write(yaml.dump(blank_database_config.model_dump(mode="json")))

    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()
    _ = runner.invoke(cli, ["db", "--config-file", str(config_path), "upgrade", "--revision", "1a9bd994df1b"])

    return config_path


def test_all_migrations(blank_initial_database_config_path):
    """Database migrations should work all the way from the oldest supported to the latest version."""
    # add some test data
    config = DbConfig.from_path(blank_initial_database_config_path)
    tan_g = "a2b6c3d9e8f7123456789abcdef0123456789abcdef0123456789abcdef01234"
    pseudonym = "CASE12345"
    submission_id = "123456789_2024-11-08_d0f805c5"
    with sqlite3.connect(config.db.database_url[len("sqlite:///") :]) as connection:
        connection.execute(
            "INSERT INTO submissions(tan_g, pseudonym, id) VALUES(:tan_g, :pseudonym, :id)",
            {"tan_g": tan_g, "pseudonym": pseudonym, "id": submission_id},
        )

    # ensure db command raises appropriate error before migration
    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()
    args_common = ["db", "--config-file", blank_initial_database_config_path]
    result_premature_list = runner.invoke(cli, [*args_common, "list"])
    assert result_premature_list.exit_code != 0
    assert "Database not at latest schema" in result_premature_list.stderr

    # run the migration
    result_upgrade = runner.invoke(cli, [*args_common, "upgrade"])
    assert result_upgrade.exit_code == 0, result_upgrade.stderr

    # check the test data
    result_show = runner.invoke(cli, [*args_common, "submission", "show", submission_id])
    assert result_show.exit_code == 0, result_show.stderr
    # shorter than tanG and less likely to be truncated in various terminal widths
    assert pseudonym in result_show.stdout, result_show.stdout


@pytest.fixture
def blank_database_config_path(tmp_path: Path, blank_database_config: DbConfig) -> Path:
    config_path = tmp_path / "config.db.yaml"
    with open(config_path, "w") as config_file:
        config_file.write(yaml.dump(blank_database_config.model_dump(mode="json")))

    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()
    _ = runner.invoke(cli, ["db", "--config-file", str(config_path), "init"])

    return config_path


def test_populate(blank_database_config_path: Path):
    args_common = ["db", "--config-file", blank_database_config_path]
    metadata = GrzSubmissionMetadata.model_validate_json(
        (importlib.resources.files(test_resources) / "metadata.json").read_text()
    )

    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()
    result_add = runner.invoke(cli, [*args_common, "submission", "add", metadata.submission_id])
    assert result_add.exit_code == 0, result_add.stderr

    with importlib.resources.as_file(importlib.resources.files(test_resources) / "metadata.json") as metadata_path:
        result_populate = runner.invoke(
            cli, [*args_common, "submission", "populate", metadata.submission_id, str(metadata_path), "--no-confirm"]
        )
    assert result_populate.exit_code == 0, result_populate.stderr

    result_show = runner.invoke(cli, [*args_common, "submission", "show", metadata.submission_id])
    assert result_show.exit_code == 0, result_show.stderr
    # shorter than tanG and less likely to be truncated in various terminal widths
    assert metadata.submission.local_case_id in result_show.stdout, result_show.stdout


def test_populate_redacted(tmp_path: Path, blank_database_config_path: Path):
    args_common = ["db", "--config-file", blank_database_config_path]
    metadata = GrzSubmissionMetadata.model_validate_json(
        (importlib.resources.files(test_resources) / "metadata.json").read_text()
    )

    # compute submission ID _before_ tanG is redacted (changing the property return value)
    submission_id = metadata.submission_id

    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()
    result_add = runner.invoke(cli, [*args_common, "submission", "add", submission_id])
    assert result_add.exit_code == 0, result_add.stderr

    # redact the tanG
    metadata.submission.tan_g = REDACTED_TAN
    metadata_path = tmp_path / "metadata.json"
    with open(metadata_path, "w") as metadata_file:
        metadata_file.write(metadata.model_dump_json(by_alias=True))

    with pytest.raises(ValueError, match="Refusing to populate a seemingly-redacted TAN"):
        _ = runner.invoke(
            cli,
            [*args_common, "submission", "populate", submission_id, str(metadata_path), "--no-confirm"],
            catch_exceptions=False,
        )


def test_update_error_confirm(blank_database_config_path: Path):
    """Database should confirm before updating a submission from an Error state."""
    args_common = ["db", "--config-file", blank_database_config_path]
    metadata = GrzSubmissionMetadata.model_validate_json(
        (importlib.resources.files(test_resources) / "metadata.json").read_text()
    )

    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()
    result_add = runner.invoke(cli, [*args_common, "submission", "add", metadata.submission_id])
    assert result_add.exit_code == 0, result_add.stderr

    result_update1 = runner.invoke(cli, [*args_common, "submission", "update", metadata.submission_id, "Error"])
    assert result_update1.exit_code == 0, result_update1.output

    result_update2 = runner.invoke(cli, [*args_common, "submission", "update", metadata.submission_id, "Validated"])
    assert result_update2.exit_code != 0, result_update2.output
    assert (
        "Submission is currently in an 'Error' state. Are you sure you want to set it to 'Validated'?"
        in result_update2.output
    )

    result_update3 = runner.invoke(
        cli, [*args_common, "submission", "update", "--ignore-error-state", metadata.submission_id, "Validated"]
    )
    assert result_update3.exit_code == 0, result_update3.output
