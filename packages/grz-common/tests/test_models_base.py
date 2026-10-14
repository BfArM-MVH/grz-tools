"""to_yaml writes secrets in plain text, so that from_path loads the file back unchanged."""

from pathlib import Path

import pytest
from grz_common.models.base import IgnoringBaseModel, IgnoringBaseSettings
from pydantic import AnyHttpUrl, SecretStr

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


class _ModelWithUrl(IgnoringBaseModel):
    url: AnyHttpUrl


@pytest.mark.filterwarnings("error")
@pytest.mark.parametrize("reveal_secrets", [False, True])
def test_json_dump_of_a_required_url_does_not_warn(reveal_secrets: bool):
    """Revealing the secrets must leave the serialization of the other fields to pydantic."""
    model = _ModelWithUrl(url="https://example.org")

    assert model.model_dump(mode="json", context={"reveal_secrets": reveal_secrets}) == {"url": "https://example.org/"}
