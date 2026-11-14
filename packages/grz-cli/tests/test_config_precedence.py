"""A ``grz_`` environment variable overrides the config file of grz-cli, as it does in grzctl."""

from grz_cli.models.config import UploadConfig, ValidateConfig
from grz_common.models.base import get_secret_value


def test_an_env_var_overrides_the_s3_credentials_of_the_file(monkeypatch):
    monkeypatch.setenv("grz_s3__access_key", "env-key")
    monkeypatch.setenv("grz_s3__secret", "env-secret")
    file_s3 = {"bucket": "test-bucket", "access_key": "file-key", "secret": "file-secret"}

    s3 = UploadConfig.model_validate({"s3": file_s3}).s3

    assert s3.access_key == "env-key"
    assert get_secret_value(s3.secret) == "env-secret"
    assert s3.bucket == "test-bucket"


def test_an_env_var_overrides_the_identifiers_of_the_file(monkeypatch):
    monkeypatch.setenv("grz_identifiers__le", "999999999")

    identifiers = ValidateConfig.model_validate({"identifiers": {"grz": "GRZX00000", "le": "123456789"}}).identifiers

    assert identifiers.le == "999999999"
    assert identifiers.grz == "GRZX00000"
