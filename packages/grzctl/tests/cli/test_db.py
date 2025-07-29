"""
Tests for grzctl db subcommand
"""

import sqlite3
from textwrap import dedent

import click.testing
import cryptography.hazmat.primitives.serialization as cryptser
import grzctl.cli
import pytest
from cryptography.hazmat.primitives.asymmetric.ed25519 import Ed25519PrivateKey
from grz_db.models.submission import OutdatedDatabaseSchemaError
from grzctl.models.config import DbConfig

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
def blank_initial_database_config_path(tmp_path):
    database_path = tmp_path / "submission.db.sqlite"
    conn = sqlite3.connect(database_path)
    conn.executescript(_INITIAL_SCHEMA)
    conn.close()

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

    config_path = tmp_path / "config.db.yaml"
    with open(config_path, "w") as config_file:
        config_file.write(
            dedent(f"""\
              db:
                database_url: "sqlite:///{database_path.resolve()}"
                author:
                  name: "alice"
                  private_key: "{private_key_path.resolve()}"
                known_public_keys: "{public_key_path.resolve()}"
        """)
        )

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
