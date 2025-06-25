import datetime
import logging
from typing import Any, ClassVar, Generic, TypeVar

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

    model_config = ConfigDict(
        json_encoders={datetime.datetime: serialize_datetime_to_iso_z},
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
    1. Define `payload_model_class: ClassVar[type[P]]`.
    2. Have an instance attribute `signature: str`.
    """

    signature: str
    payload_model_class: ClassVar[type[P]]

    def __init_subclass__(cls, **kwargs: Any) -> None:  # noqa: D105
        super().__init_subclass__(**kwargs)
        if not hasattr(cls, "payload_model_class"):
            raise TypeError(f"Class {cls.__name__} lacks 'payload_model_class' attribute required by VerifiableLog.")
        if not (isinstance(cls.payload_model_class, type) and issubclass(cls.payload_model_class, BaseSignablePayload)):
            raise TypeError(
                f"'payload_model_class' in {cls.__name__} must be a class and a subclass of BaseSignedPayload. "
                f"Got: {cls.payload_model_class}"
            )

    def verify(self, public_key: Ed25519PublicKey) -> bool:
        """Verify the signature of this log entry."""
        if not hasattr(self, "signature") or not isinstance(self.signature, str) or not self.signature:
            log.warning(f"Missing/invalid signature for {self.__class__.__name__} (id: {getattr(self, 'id', 'N/A')}).")
            return False

        signature_bytes = bytes.fromhex(self.signature)
        data_for_payload = self.model_dump(by_alias=True, exclude={"signature", "payload_model_class"})
        payload_to_verify = self.payload_model_class(**data_for_payload)
        bytes_to_verify = payload_to_verify.to_bytes()

        try:
            public_key.verify(signature_bytes, bytes_to_verify)
        except cryptography.exceptions.InvalidSignature:
            return False
        except:
            raise
        return True
