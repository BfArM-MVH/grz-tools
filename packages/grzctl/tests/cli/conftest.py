import importlib.resources
from pathlib import Path

import click.testing
import cryptography.hazmat.primitives.serialization as cryptser
import grzctl.cli
import pytest
import yaml
from cryptography.hazmat.primitives.asymmetric.ed25519 import Ed25519PrivateKey
from grz_pydantic_models_testing.example_metadata import grzctl as grzctl_metadata
from grzctl.models.config import DbConfig

#: The revision the schema-upgrade tests start from.
INITIAL_REVISION = "1a9bd994df1b"


@pytest.fixture
def test_metadata_path():
    return importlib.resources.files(grzctl_metadata).joinpath("metadata.json")


def _database_config(tmp_path: Path, database_url: str) -> DbConfig:
    """Build a config for *database_url*, signed by a throwaway key pair written under *tmp_path*."""
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
            "database_url": database_url,
            "author": {
                "name": "alice",
                "private_key_path": str(private_key_path.resolve()),
                "private_key_passphrase": "",
            },
            "known_public_keys": str(public_key_path.resolve()),
        }
    )


def _write_config(tmp_path: Path, config: DbConfig) -> Path:
    config_path = tmp_path / "config.db.yaml"
    with open(config_path, "w") as config_file:
        config_file.write(yaml.dump(config.model_dump(mode="json")))
    return config_path


@pytest.fixture
def blank_database_config(tmp_path: Path, migrated_db_connection: str) -> DbConfig:
    """Config for a database already on the latest schema, one per supported backend."""
    return _database_config(tmp_path, migrated_db_connection)


@pytest.fixture
def blank_database_config_path(tmp_path: Path, blank_database_config: DbConfig) -> Path:
    """Path to a config file for a database on the latest schema."""
    return _write_config(tmp_path, blank_database_config)


@pytest.fixture
def empty_database_config_path(tmp_path: Path, db_test_connection: str) -> Path:
    """Path to a config file for a database with no schema at all.

    For tests that drive ``db init`` or ``db upgrade`` themselves.
    """
    return _write_config(tmp_path, _database_config(tmp_path, db_test_connection))


@pytest.fixture
def blank_initial_database_config_path(empty_database_config_path: Path) -> Path:
    """Path to a config file for a database at :data:`INITIAL_REVISION`, ready to be upgraded."""
    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()
    result = runner.invoke(
        cli, ["db", "--config-file", str(empty_database_config_path), "upgrade", "--revision", INITIAL_REVISION]
    )
    assert result.exit_code == 0, result.stderr

    return empty_database_config_path
