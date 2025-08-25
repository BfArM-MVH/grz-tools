from datetime import UTC, datetime, timedelta
from typing import Annotated

import jwt
from fastapi import Depends, HTTPException, status
from fastapi.security import OAuth2PasswordBearer
from jwt import PyJWTError
from passlib.context import CryptContext
from pydantic import BaseModel

from .config import GatekeeperConfig, get_gatekeeper_config

pwd_context = CryptContext(schemes=["bcrypt"], deprecated="auto")
oauth2_scheme = OAuth2PasswordBearer(tokenUrl="/v1/token")


class Token(BaseModel):
    access_token: str
    token_type: str


class TokenData(BaseModel):
    username: str | None = None


class AuthUser(BaseModel):
    username: str
    disabled: bool | None = None


class UserInDB(AuthUser):
    hashed_password: str


def get_user_from_config(config: GatekeeperConfig, username: str) -> UserInDB | None:
    """Fetches user details from the loaded configuration."""
    if username in config.auth.users:
        user_config = config.auth.users[username]
        return UserInDB(username=username, hashed_password=user_config.hashed_password, disabled=user_config.disabled)
    return None


def verify_password(plain_password: str, hashed_password: str) -> bool:
    return pwd_context.verify(plain_password, hashed_password)


def authenticate_user(config: GatekeeperConfig, username: str, password: str) -> AuthUser | None:
    user = get_user_from_config(config, username)
    if not user:
        return None
    if not verify_password(password, user.hashed_password):
        return None
    return AuthUser(username=user.username, disabled=user.disabled)


def create_access_token(data: dict, config: GatekeeperConfig) -> str:
    to_encode = data.copy()
    expire = datetime.now(UTC) + timedelta(minutes=config.auth.access_token_expire_minutes)
    to_encode.update({"exp": expire})
    encoded_jwt = jwt.encode(to_encode, config.auth.secret_key, algorithm=config.auth.algorithm)
    return encoded_jwt


async def get_current_user(
    token: Annotated[str, Depends(oauth2_scheme)], config: GatekeeperConfig = Depends(get_gatekeeper_config)
) -> AuthUser:
    credentials_exception = HTTPException(
        status_code=status.HTTP_401_UNAUTHORIZED,
        detail="Could not validate credentials",
        headers={"WWW-Authenticate": "Bearer"},
    )
    try:
        payload = jwt.decode(token, config.auth.secret_key, algorithms=[config.auth.algorithm])
        username: str | None = payload.get("sub")
        if username is None:
            raise credentials_exception
        token_data = TokenData(username=username)
    except PyJWTError as e:
        raise credentials_exception from e

    user_config = get_user_from_config(config, token_data.username)
    if user_config is None:
        raise credentials_exception
    return AuthUser(username=user_config.username, disabled=user_config.disabled)


async def get_current_active_user(current_user: Annotated[AuthUser, Depends(get_current_user)]) -> AuthUser:
    if current_user.disabled:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="Inactive user")
    return current_user
