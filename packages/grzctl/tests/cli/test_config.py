"""GrzctlConfig checks the configuration it loads, and merges environment variables into it."""

import json

import pytest
from grz_common.models.base import get_secret_value
from grzctl.models.config import GrzctlConfig
from pydantic import ValidationError

LE_ID = "260914050"
BUCKET_NAME = "grz-inbox-test"
INBOX = {
    "endpoint_url": "https://s3.amazonaws.com",
    "access_key": "testing",
    "secret": "testing",
}


@pytest.fixture
def configuration(offline_config: GrzctlConfig) -> dict:
    """The offline config as a dict, whose only inbox is ``BUCKET_NAME`` of submitter ``LE_ID``."""
    configuration = offline_config.model_dump(mode="json", exclude_none=True)
    inbox = {**INBOX, "private_key_path": configuration["archives"]["signing_key_path"]}
    configuration["leistungserbringer"] = {LE_ID: {"inbox_buckets": {BUCKET_NAME: inbox}}}
    return configuration


def test_pydantic_nested_env_var_merging(monkeypatch, configuration: dict):
    env_var_name = f"GRZ_LEISTUNGSERBRINGER__{LE_ID}__INBOX_BUCKETS__{BUCKET_NAME}__PRIVATE_KEY_PASSPHRASE"
    monkeypatch.setenv(env_var_name.upper(), "dotenv-secret-passphrase")

    config = GrzctlConfig.from_configuration(configuration)

    entry = config.leistungserbringer[LE_ID]
    assert get_secret_value(entry.inbox_buckets[BUCKET_NAME].private_key_passphrase) == "dotenv-secret-passphrase"


def test_pydantic_json_env_var_merging(monkeypatch, configuration: dict):
    inbox_override = {**INBOX, "private_key_passphrase": "json-secret-passphrase"}
    monkeypatch.setenv("GRZ_LEISTUNGSERBRINGER", json.dumps({LE_ID: {"inbox_buckets": {BUCKET_NAME: inbox_override}}}))

    config = GrzctlConfig.from_configuration(configuration)

    entry = config.leistungserbringer[LE_ID]
    assert get_secret_value(entry.inbox_buckets[BUCKET_NAME].private_key_passphrase) == "json-secret-passphrase"


def test_archive_public_key_can_come_from_an_env_var(monkeypatch, configuration: dict, crypt4gh_public_key: str):
    """An operator may put the archive's public key into an environment variable inline,
    instead of writing it to a file that ``public_key_path`` then points at.
    """
    del configuration["archives"]["consented"]["public_key_path"]
    monkeypatch.setenv("GRZ_ARCHIVES__CONSENTED__PUBLIC_KEY", crypt4gh_public_key)

    config = GrzctlConfig.from_configuration(configuration)

    assert config.archives.consented.public_key == crypt4gh_public_key
    assert config.archives.consented.public_key_path is None


@pytest.mark.parametrize("field", ["authorization_url", "client_id", "client_secret", "api_base_url"])
def test_a_config_without_a_pruefbericht_field_fails(configuration: dict, field: str):
    """Every grzctl command loads the whole config.
    So a missing field stops even the commands that submit no Prüfbericht.
    """
    del configuration["pruefbericht"][field]

    with pytest.raises(ValidationError, match=rf"pruefbericht\.{field}\n  Field required"):
        GrzctlConfig.from_configuration(configuration)
