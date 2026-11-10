"""to_yaml writes secrets in plain text, so that from_path loads the file back unchanged."""

from pathlib import Path

import pytest
from grz_common.models.base import IgnoringBaseModel, IgnoringBaseSettings
from pydantic import SecretStr

SECRET = "s3-secret"


class _Model(IgnoringBaseModel):
    secret: SecretStr


class _Settings(IgnoringBaseSettings):
    secret: SecretStr


@pytest.mark.parametrize("model_class", [_Model, _Settings])
def test_to_yaml_roundtrips_secrets(tmp_path: Path, model_class: type[_Model | _Settings]):
    config_path = tmp_path / "config.yaml"
    with open(config_path, "w") as fd:
        model_class(secret=SecretStr(SECRET)).to_yaml(fd)

    assert model_class.from_path(config_path).secret.get_secret_value() == SECRET
