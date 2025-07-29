"""
Tests for grzctl db subcommand
"""

import hashlib
import importlib.resources
import json
import sqlite3
from pathlib import Path

import click.testing
import cryptography.hazmat.primitives.serialization as cryptser
import grzctl.cli
import pytest
import yaml
from cryptography.hazmat.primitives.asymmetric.ed25519 import Ed25519PrivateKey
from grz_db.models.submission import OutdatedDatabaseSchemaError
from grz_pydantic_models.submission.metadata import GrzSubmissionMetadata
from grzctl.models.config import DbConfig

from .. import resources as test_resources

_INITIAL_SCHEMA = """\
CREATE TABLE submissions (
 tan_g VARCHAR,
 pseudonym VARCHAR,
 id VARCHAR NOT NULL,
 PRIMARY KEY (id)
);
CREATE UNIQUE INDEX ix_submissions_tan_g ON submissions (tan_g);
CREATE INDEX ix_submissions_id ON submissions (id);
CREATE INDEX ix_submissions_pseudonym ON submissions (pseudonym);
CREATE TABLE submission_states (
 state VARCHAR(11) NOT NULL,
 data JSON,
 timestamp DATETIME NOT NULL,
 id INTEGER NOT NULL,
 submission_id VARCHAR NOT NULL,
 author_name VARCHAR NOT NULL,
 signature VARCHAR NOT NULL,
 PRIMARY KEY (id),
 FOREIGN KEY(submission_id) REFERENCES submissions (id)
);
CREATE INDEX ix_submission_states_author_name ON submission_states (author_name);
CREATE INDEX ix_submission_states_submission_id ON submission_states (submission_id);
CREATE INDEX ix_submission_states_id ON submission_states (id);
CREATE TABLE submission_change_requests (
 change VARCHAR(8) NOT NULL,
 data JSON,
 timestamp DATETIME NOT NULL,
 id INTEGER NOT NULL,
 submission_id VARCHAR NOT NULL,
 author_name VARCHAR NOT NULL,
 signature VARCHAR NOT NULL,
 PRIMARY KEY (id),
 FOREIGN KEY(submission_id) REFERENCES submissions (id)
);
CREATE INDEX ix_submission_change_requests_author_name ON submission_change_requests (author_name);
CREATE INDEX ix_submission_change_requests_id ON submission_change_requests (id);
CREATE INDEX ix_submission_change_requests_submission_id ON submission_change_requests (submission_id);
"""


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
            "author": {"name": "alice", "private_key_path": str(private_key_path.resolve())},
            "known_public_keys": str(public_key_path.resolve()),
        }
    )


@pytest.fixture
def blank_initial_database_config_path(tmp_path: Path, blank_database_config: DbConfig) -> Path:
    conn = sqlite3.connect(blank_database_config.db.database_url[len("sqlite:///") :])
    conn.executescript(_INITIAL_SCHEMA)
    conn.close()

    config_path = tmp_path / "config.db.yaml"
    with open(config_path, "w") as config_file:
        config_file.write(yaml.dump(blank_database_config.model_dump(mode="json")))

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
    with pytest.raises(OutdatedDatabaseSchemaError):
        _ = runner.invoke(cli, [*args_common, "list"], catch_exceptions=False)

    # run the migration
    result_upgrade = runner.invoke(cli, [*args_common, "upgrade"])
    assert result_upgrade.exit_code == 0, result_upgrade.stderr

    # check the test data
    result_show = runner.invoke(cli, [*args_common, "submission", "show", submission_id])
    assert result_show.exit_code == 0, result_show.stderr
    assert f"tanG: {tan_g}" in result_show.stdout


@pytest.fixture
def blank_database_config_path(tmp_path: Path, blank_database_config: DbConfig) -> Path:
    config_path = tmp_path / "config.db.yaml"
    with open(config_path, "w") as config_file:
        config_file.write(yaml.dump(blank_database_config.model_dump(mode="json")))

    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()
    _ = runner.invoke(cli, ["db", "--config-file", str(config_path), "init"])

    return config_path


def _generate_submission_id(metadata: GrzSubmissionMetadata) -> str:
    # generate a submission ID once
    submitter_id = metadata.submission.submitter_id
    submission_date = metadata.submission.submission_date
    # use first 8 characters of SHA256 hash of transaction ID to virtually prevent collisions
    suffix = hashlib.sha256(metadata.submission.tan_g.encode("utf-8")).hexdigest()[:8]
    return f"{submitter_id}_{submission_date}_{suffix}"


def test_populate(blank_database_config_path: Path):
    args_common = ["db", "--config-file", blank_database_config_path]
    metadata = GrzSubmissionMetadata.model_validate_json(
        (importlib.resources.files(test_resources) / "metadata.json").read_text()
    )
    submission_id = _generate_submission_id(metadata)

    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()
    result_add = runner.invoke(cli, [*args_common, "submission", "add", submission_id])
    assert result_add.exit_code == 0, result_add.stderr

    with importlib.resources.as_file(importlib.resources.files(test_resources) / "metadata.json") as metadata_path:
        result_populate = runner.invoke(
            cli, [*args_common, "submission", "populate", submission_id, str(metadata_path), "--accept-changes"]
        )
    assert result_populate.exit_code == 0, result_populate.stderr

    result_show = runner.invoke(cli, [*args_common, "submission", "show", submission_id])
    assert result_show.exit_code == 0, result_show.stderr
    assert f"tanG: {metadata.submission.tan_g}" in result_show.stdout


def test_populate_redacted(tmp_path: Path, blank_database_config_path: Path):
    args_common = ["db", "--config-file", blank_database_config_path]
    metadata = GrzSubmissionMetadata.model_validate_json(
        (importlib.resources.files(test_resources) / "metadata.json").read_text()
    )
    submission_id = _generate_submission_id(metadata)

    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()
    result_add = runner.invoke(cli, [*args_common, "submission", "add", submission_id])
    assert result_add.exit_code == 0, result_add.stderr

    # redact the tanG
    metadata.submission.tan_g = "0" * 64
    metadata_path = tmp_path / "metadata.json"
    with open(metadata_path, "w") as metadata_file:
        json.dump(metadata.model_dump(mode="json"), metadata_file)

    with pytest.raises(ValueError):
        _ = runner.invoke(
            cli,
            [*args_common, "submission", "populate", submission_id, str(metadata_path), "--accept-changes"],
            catch_exceptions=False,
        )
