from typing import Annotated

from grz_common.models.base import IgnoringBaseModel
from pydantic import AnyHttpUrl, SecretStr, UrlConstraints


class PruefberichtModel(IgnoringBaseModel):
    authorization_url: Annotated[AnyHttpUrl, UrlConstraints(allowed_schemes=["https"], host_required=True)] | None = (
        None
    )
    """
    URL from which to request a new Prüfbericht submission token
    """

    client_id: str | None = None
    """
    Client ID used to obtain new Prüfbericht submission tokens
    """

    client_secret: SecretStr | None = None
    """
    Client secret used to obtain new Prüfbericht submission tokens
    """

    api_base_url: Annotated[AnyHttpUrl, UrlConstraints(allowed_schemes=["https"], host_required=True)] | None = None
    """
    Base URL to BfArM Submission (Prüfbericht) API
    """
