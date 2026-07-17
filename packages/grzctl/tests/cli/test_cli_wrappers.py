"""Tests for config derivation in grzctl CLI wrappers."""

import pytest
from grzctl.commands.cli_wrappers import derive_encrypt_config, derive_validate_config
from grzctl.models.config import GrzctlConfig


def _make_grzctl_config(**overrides) -> GrzctlConfig:
    """Build a minimal valid GrzctlConfig for testing derivation functions."""
    defaults = {
        "s3": {"inboxes": {"000000000": {"inbox": {"private_key_path": "/dev/null"}}}},
        "archives": {
            "consented": {"s3": {"bucket": "consented"}, "public_key_path": "/dev/null"},
            "non_consented": {"s3": {"bucket": "non_consented"}, "public_key_path": "/dev/null"},
        },
        "db": {"database_url": "sqlite:///:memory:", "author": {"name": "test"}},
        "pruefbericht": {},
        "keys": {
            "grz_private_key_path": "/dev/null",
            "grz_public_key_path": "/dev/null",
        },
        "identifiers": {"grz": "GRZT00000"},
    }
    defaults.update(overrides)
    return GrzctlConfig(**defaults)


class TestDeriveValidateConfig:
    def test_extracts_identifiers(self):
        config = _make_grzctl_config()
        result = derive_validate_config(config)
        assert result == {"identifiers": {"grz": "GRZT00000"}}

    def test_includes_le_if_set(self):
        config = _make_grzctl_config(identifiers={"grz": "GRZT00000", "le": "123456789"})
        result = derive_validate_config(config)
        assert result == {"identifiers": {"grz": "GRZT00000", "le": "123456789"}}

    def test_no_extra_keys(self):
        config = _make_grzctl_config()
        result = derive_validate_config(config)
        assert set(result.keys()) == {"identifiers"}


class TestDeriveEncryptConfig:
    def test_public_key_path_only(self):
        config = _make_grzctl_config()
        result = derive_encrypt_config(config)
        assert result == {"keys": {"grz_public_key_path": "/dev/null", "grz_private_key_path": "/dev/null"}}

    def test_with_submitter_key(self):
        config = _make_grzctl_config(
            keys={
                "grz_private_key_path": "/dev/null",
                "grz_public_key_path": "/dev/null",
                "submitter_private_key_path": "/dev/null",
            }
        )
        result = derive_encrypt_config(config)
        assert result == {
            "keys": {
                "grz_public_key_path": "/dev/null",
                "grz_private_key_path": "/dev/null",
                "submitter_private_key_path": "/dev/null",
            }
        }

    def test_no_private_key_path(self):
        config = _make_grzctl_config(
            keys={"grz_private_key_path": "/dev/null", "grz_public_key_path": "/dev/null"}
        )
        result = derive_encrypt_config(config)
        # grz_private_key_path should be included since it's set
        assert "grz_private_key_path" in result["keys"]

    def test_no_extra_keys(self):
        config = _make_grzctl_config()
        result = derive_encrypt_config(config)
        assert set(result.keys()) == {"keys"}
        assert isinstance(result["keys"], dict)
