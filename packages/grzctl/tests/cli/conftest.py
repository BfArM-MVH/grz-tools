from pathlib import Path

import click.testing
import cryptography.hazmat.primitives.serialization as cryptser
import grzctl.cli
import pytest
import yaml
from cryptography.hazmat.primitives.asymmetric.ed25519 import Ed25519PrivateKey
from grzctl.models.config import DbConfig


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


@pytest.fixture
def blank_database_config_path(tmp_path: Path, blank_database_config: DbConfig) -> Path:
    config_path = tmp_path / "config.db.yaml"
    with open(config_path, "w") as config_file:
        config_file.write(yaml.dump(blank_database_config.model_dump(mode="json")))

    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()
    _ = runner.invoke(cli, ["db", "--config-file", str(config_path), "init"])

    return config_path
