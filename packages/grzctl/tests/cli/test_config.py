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
    inbox = {**INBOX, "private_key_path": configuration["db"]["author"]["private_key_path"]}
    configuration["leistungserbringer"] = {LE_ID: {"inbox_buckets": {BUCKET_NAME: inbox}}}
    return configuration


def test_pydantic_nested_env_var_merging(monkeypatch, configuration: dict):
    env_var_name = f"GRZ_LEISTUNGSERBRINGER__{LE_ID}__INBOX_BUCKETS__{BUCKET_NAME}__PRIVATE_KEY_PASSPHRASE"
    monkeypatch.setenv(env_var_name.upper(), "dotenv-secret-passphrase")

    config = GrzctlConfig.from_configuration(configuration)

    entry = config.leistungserbringer[LE_ID]
    assert get_secret_value(entry.inbox_buckets[BUCKET_NAME].private_key_passphrase) == "dotenv-secret-passphrase"


def test_an_env_var_reaches_an_inbox_whose_name_has_upper_case_letters(monkeypatch, configuration: dict):
    """pydantic-settings reads the names of environment variables in lowercase.
    The inbox in a variable's name still matches the inbox of the file whose name differs only in case.
    """
    inboxes = configuration["leistungserbringer"][LE_ID]["inbox_buckets"]
    inboxes["Main"] = inboxes.pop(BUCKET_NAME)
    monkeypatch.setenv(f"GRZ_LEISTUNGSERBRINGER__{LE_ID}__INBOX_BUCKETS__MAIN__SECRET", "env-secret")

    config = GrzctlConfig.from_configuration(configuration)

    inbox_buckets = config.leistungserbringer[LE_ID].inbox_buckets
    assert list(inbox_buckets) == ["Main"]
    assert get_secret_value(inbox_buckets["Main"].secret) == "env-secret"


def test_inbox_names_must_differ_in_more_than_case(configuration: dict):
    """An environment variable names an inbox in lowercase, so it could not tell these inboxes apart."""
    inboxes = configuration["leistungserbringer"][LE_ID]["inbox_buckets"]
    inboxes["Main"] = inboxes["main"] = inboxes.pop(BUCKET_NAME)

    with pytest.raises(ValidationError, match=r"Inbox names must differ in more than case: 'Main', 'main'"):
        GrzctlConfig.from_configuration(configuration)


def test_two_env_vars_that_spell_one_inbox_differently_fail(monkeypatch, configuration: dict):
    """A JSON variable keeps the case of the inbox name, while a nested variable gives it in lowercase.
    The two must not override each other silently.
    """
    inboxes = configuration["leistungserbringer"][LE_ID]["inbox_buckets"]
    inboxes["Main"] = inboxes.pop(BUCKET_NAME)
    monkeypatch.setenv(f"GRZ_LEISTUNGSERBRINGER__{LE_ID}__INBOX_BUCKETS", json.dumps({"Main": {"secret": "json"}}))
    monkeypatch.setenv(f"GRZ_LEISTUNGSERBRINGER__{LE_ID}__INBOX_BUCKETS__MAIN__ACCESS_KEY", "nested")

    with pytest.raises(ValidationError, match=r"Inbox names must differ in more than case: 'Main', 'main'"):
        GrzctlConfig.from_configuration(configuration)


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


@pytest.mark.parametrize("env_var_name", ["GRZ_PRIVATE_KEY", "GRZ_PRIVATE_KEY_PASSPHRASE"])
def test_db_author_ignores_env_vars_without_the_db_author_path(monkeypatch, configuration: dict, env_var_name: str):
    """An environment variable sets a field of ``db.author`` only by its full path,
    such as ``GRZ_DB__AUTHOR__PRIVATE_KEY_PASSPHRASE``.
    """
    configuration["db"]["author"].pop("private_key_passphrase", None)
    monkeypatch.setenv(env_var_name, "stray-value")

    config = GrzctlConfig.from_configuration(configuration)

    assert config.db.author.private_key is None
    assert config.db.author.private_key_passphrase is None


def test_inbox_target_defaults_the_bucket_to_the_inbox_name(configuration: dict):
    """Without an explicit ``bucket:``, the S3 bucket of an inbox is its name."""
    config = GrzctlConfig.from_configuration(configuration)

    target = config.inbox_target(LE_ID, BUCKET_NAME)

    assert target.s3.bucket == BUCKET_NAME


def test_inbox_target_honors_an_explicit_bucket_override(configuration: dict):
    """An explicit ``bucket:`` names an S3 bucket that differs from the inbox name."""
    inbox = {
        **INBOX,
        "private_key_path": configuration["db"]["author"]["private_key_path"],
        "bucket": "grz-incoming-prod",
    }
    configuration["leistungserbringer"] = {LE_ID: {"inbox_buckets": {"inbox-external": inbox}}}

    config = GrzctlConfig.from_configuration(configuration)

    assert config.inbox_target(LE_ID, "inbox-external").s3.bucket == "grz-incoming-prod"


@pytest.mark.parametrize("field", ["authorization_url", "client_id", "client_secret", "api_base_url"])
def test_a_config_without_a_pruefbericht_field_fails(configuration: dict, field: str):
    """Every grzctl command loads the whole config.
    So a missing field stops even the commands that submit no Prüfbericht.
    """
    del configuration["pruefbericht"][field]

    with pytest.raises(ValidationError, match=rf"pruefbericht\.{field}\n  Field required"):
        GrzctlConfig.from_configuration(configuration)
