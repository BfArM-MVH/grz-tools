import datetime
import logging
from typing import Any, Generic, TypeVar

import cryptography
from cryptography.hazmat.primitives.asymmetric.ed25519 import Ed25519PrivateKey, Ed25519PublicKey
from pydantic import ConfigDict
from sqlmodel import SQLModel

from ..common import serialize_datetime_to_iso_z

log = logging.getLogger(__name__)


class BaseSignablePayload(SQLModel):
    """
    Base class for SQLModel based payloads
    that can be signed and can be converted to bytes for verification.
    Provides a default `to_bytes` method using pydantic's JSON serialization.
    Provides a default `sign` method using the private key of the author.
    """

    model_config = ConfigDict(  # type: ignore
        populate_by_name=True,
    )

    def to_bytes(self) -> bytes:
        """
        Default serialization: JSON string encoded to UTF-8.
        """
        payload_json = self.model_dump_json(by_alias=True)
        return payload_json.encode("utf8")

    def sign(self, private_key: Ed25519PrivateKey) -> bytes:
        """Sign this payload using the given private key."""
        bytes_to_sign = self.to_bytes()
        signature = private_key.sign(bytes_to_sign)
        public_key_of_private = private_key.public_key()
        public_key_of_private.verify(signature, bytes_to_sign)
        return signature


P = TypeVar("P", bound=BaseSignablePayload)


class VerifiableLog(Generic[P]):
    """
    Mixin class for SQLModels that store a signature and can be verified.
    Subclasses MUST:
    1. Define `_payload_model_class: type[P]`.
    2. Have an instance attribute `signature: str`.
    """

    signature: str
    _payload_model_class: type[P]

    def __init_subclass__(cls, **kwargs: Any) -> None:  # noqa: D105
        super().__init_subclass__(**kwargs)
        if not hasattr(cls, "_payload_model_class"):
            raise TypeError(f"Class {cls.__name__} lacks '_payload_model_class' attribute required by VerifiableLog.")
        payload_class_or_field = getattr(cls, "_payload_model_class")

        if hasattr(payload_class_or_field, "default"):
            payload_class = payload_class_or_field.default
        else:
            payload_class = payload_class_or_field

        if not (isinstance(payload_class, type) and issubclass(payload_class, BaseSignablePayload)):
            raise TypeError(
                f"'_payload_model_class' in {cls.__name__} must be a class and a subclass of BaseSignedPayload. "
                f"Got: {payload_class_or_field}"
            )

    def verify(self, public_key: Ed25519PublicKey) -> bool:
        """Verify the signature of this log entry."""
        if not hasattr(self, "signature") or not isinstance(self.signature, str) or not self.signature:
            log.warning(f"Missing/invalid signature for {self.__class__.__name__} (id: {getattr(self, 'id', 'N/A')}).")
            return False

        signature_bytes = bytes.fromhex(self.signature)

        payload_fields = self._payload_model_class.model_fields

        data_for_payload = {
            field_name: getattr(self, field_name) for field_name in payload_fields if hasattr(self, field_name)
        }

        payload_to_verify = self._payload_model_class(**data_for_payload)
        bytes_to_verify = payload_to_verify.to_bytes()

        try:
            public_key.verify(signature_bytes, bytes_to_verify)
        except cryptography.exceptions.InvalidSignature:
            return False
        except:
            log.error("An unexpected error occurred during signature verification:", exc_info=True)
            raise
        return True
