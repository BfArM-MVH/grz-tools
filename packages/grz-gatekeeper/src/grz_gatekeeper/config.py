import logging
from functools import lru_cache

from fastapi import HTTPException, status
from grz_common.models.base import IgnoringBaseSettings
from grz_common.models.s3 import S3ConfigModel
from grzctl.models.config import DbConfig
from pydantic import Field

log = logging.getLogger(__name__)

CONFIG_FILE_PATH: str | None = None


class AuthUserConfig(IgnoringBaseSettings):
    hashed_password: str
    disabled: bool = False


class AuthConfig(IgnoringBaseSettings):
    secret_key: str
    algorithm: str
    access_token_expire_minutes: int
    users: dict[str, AuthUserConfig] = Field(default_factory=dict)


class GatekeeperConfig(DbConfig, S3ConfigModel):
    auth: AuthConfig
    max_active_sessions: int = Field(100, ge=1)


@lru_cache
def get_gatekeeper_config() -> GatekeeperConfig:
    """Loads, validates, and caches the application configuration from a file."""
    if CONFIG_FILE_PATH is None:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Server configuration path not set.",
        )
    try:
        return GatekeeperConfig.from_path(CONFIG_FILE_PATH)
    except FileNotFoundError as e:
        log.error(f"Configuration file not found at: {CONFIG_FILE_PATH}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail="Server configuration is missing."
        ) from e
