"""GrzctlConfig merges environment variables into the configuration it loads."""

import json

import pytest
from grz_common.models.base import get_secret_value
from grzctl.models.config import GrzctlConfig

LE_ID = "260914050"
BUCKET_NAME = "grz-inbox-test"
INBOX = {
    "endpoint_url": "https://s3.amazonaws.com",
    "access_key": "testing",
    "secret": "testing",
    "private_key_path": "/path/to/test.sec",
}


@pytest.fixture
def configuration(offline_config: GrzctlConfig) -> dict:
    """The offline config as a dict, whose only inbox is ``BUCKET_NAME`` of submitter ``LE_ID``."""
    configuration = offline_config.model_dump(mode="json", exclude_none=True)
    configuration["leistungserbringer"] = {LE_ID: {"inbox_buckets": {BUCKET_NAME: INBOX}}}
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
