"""A ``grz_`` environment variable sets an S3 credential that the config file of grz-cli leaves out."""

from grz_cli.models.config import UploadConfig
from grz_common.models.base import get_secret_value

S3 = {"endpoint_url": "https://s3.example.org", "bucket": "test-bucket"}


def test_an_env_var_sets_the_s3_credential_that_the_file_leaves_out(monkeypatch):
    monkeypatch.setenv("grz_s3__secret", "env-secret")

    s3 = UploadConfig.model_validate({"s3": {**S3, "access_key": "file-key"}}).s3

    assert s3.access_key == "file-key"
    assert get_secret_value(s3.secret) == "env-secret"
