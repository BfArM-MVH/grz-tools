import importlib.resources
import shutil
from pathlib import Path

import click.testing
import cryptography.hazmat.primitives.serialization as cryptser
import grzctl.cli
import psycopg
import pytest
from cryptography.hazmat.primitives.asymmetric.ed25519 import Ed25519PrivateKey
from grz_pydantic_models_testing.example_metadata import grzctl as grzctl_metadata
from grzctl.models.config import GrzctlConfig


def _grzctl_archives(endpoint_url: str | None = None, public_key_path: str = "/dev/null") -> dict:

    def _s3(bucket):
        d = {"bucket": bucket, "public_key_path": public_key_path}
        if endpoint_url:
            d["endpoint_url"] = endpoint_url
        return d

    return {
        "consented": {"s3": _s3("consented"), "public_key_path": public_key_path},
        "non_consented": {"s3": _s3("non_consented"), "public_key_path": public_key_path},
    }


@pytest.fixture
def test_metadata_path():
    return importlib.resources.files(grzctl_metadata).joinpath("metadata.json")


@pytest.fixture(
    params=[
        "sqlite",
        pytest.param(
            "postgresql",
            marks=pytest.mark.skipif(condition=shutil.which("pg_config") is None, reason="postgresql not detected"),
        ),
    ]
)
def blank_database_config(request: pytest.FixtureRequest, tmp_path: Path) -> GrzctlConfig:
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

    database_url = "sqlite:///" + str((tmp_path / "submission.db.sqlite").resolve())
    if request.param == "postgresql":
        postgresql: psycopg.Connection = request.getfixturevalue("postgresql")
        database_url = f"postgresql+psycopg://{postgresql.info.user}:@{postgresql.info.host}:{postgresql.info.port}/{postgresql.info.dbname}"

    return GrzctlConfig(
        s3={
            "inboxes": {
                "000000000": {
                    "inbox": {
                        "private_key_path": str(private_key_path.resolve()),
                    }
                }
            }
        },
        archives=_grzctl_archives(
            public_key_path=str(public_key_path.resolve()),
        ),
        db={
            "database_url": database_url,
            "author": {
                "name": "alice",
                "private_key_path": str(private_key_path.resolve()),
                "private_key_passphrase": "",
            },
            "known_public_keys": str(public_key_path.resolve()),
        },
        pruefbericht={},
        keys={"grz_public_key_path": str(public_key_path.resolve())},
        identifiers={"grz": "GRZK00007"},
    )


@pytest.fixture
def blank_initial_database_config_path(tmp_path: Path, blank_database_config: GrzctlConfig) -> Path:
    config_path = tmp_path / "config.db.yaml"
    with open(config_path, "w") as config_file:
        blank_database_config.to_yaml(config_file)

    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()
    _ = runner.invoke(cli, ["db", "--config-file", str(config_path), "upgrade", "--revision", "1a9bd994df1b"])

    return config_path


@pytest.fixture
def blank_database_config_path(tmp_path: Path, blank_database_config: GrzctlConfig) -> Path:
    import yaml

    config_path = tmp_path / "config.db.yaml"

    config_dict = blank_database_config.model_dump(mode="json", exclude_none=True, exclude_unset=True, exclude_defaults=True)

    with open(config_path, "w") as config_file:
        yaml.dump(config_dict, config_file)

    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()
    _ = runner.invoke(cli, ["db", "--config-file", str(config_path), "init"])

    return config_path
