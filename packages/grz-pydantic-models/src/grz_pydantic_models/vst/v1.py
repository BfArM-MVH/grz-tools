"""
Models for the pseudonym interface of the RKI Vertrauensstelle (VST), API version 1.

The Vertrauensstelle of the Modellvorhaben Genomsequenzierung at the RKI issues Vorgangsnummern
(the ``tanG``/``tanC`` of a submission) and Pseudonyme. This module holds the request and
response bodies of its three REST endpoints (``POST``, JSON):

* ``psngen/api/v1/pseudonym/migrate`` (Phase 0): :class:`MigrationRequestData` in, HTTP 200 with
  an empty body out. A Leistungserbringer (LE) reports data with a Vorgangsnummer it made itself.
* ``psngen/api/v1/pseudonym/vn`` (Phase 1; the test server serves it as ``.../generate``):
  :class:`VnRequestData` in, :class:`VnResult` out. The VST issues the Vorgangsnummer.
* ``psn/api/v1/pseudonym/result`` (both phases): :class:`PseudonymRequestData` in,
  :class:`PseudonymResult` out. A data node (GRZ or KDK) fetches the Pseudonym for a Vorgangsnummer.

Errors of every endpoint (400, 401, 403, 404, 415, 500) carry an :class:`ErrorDetails` body.

Provenance
----------

The RKI publishes the interface as three OpenAPI documents, kept in
``grz_pydantic_models/resources/vst/``. They are authoritative; the PDF specification v1.0 and the
test handout lag behind them in places. This module started as the output of::

    datamodel-codegen --input <the three documents merged> --input-file-type openapi
        --output-model-type pydantic_v2.BaseModel --target-python-version 3.12
        --use-standard-collections --use-union-operator --field-constraints
        --use-schema-description --use-field-description

and was then edited by hand: the field docstrings are the OpenAPI descriptions verbatim, the
classes, field names on the wire, requiredness and ``extra="forbid"`` are the generated ones.
To pick up a revised document, regenerate into a scratch file and diff it against this module.

The hand edits add what the documents cannot express:

* Field values are validated against the GRZ metadata rules (KVNR, TIDs, hexadecimal
  Vorgangsnummer) instead of being arbitrary Base64 strings.
* ``PSN`` and ``Signature`` are ``bytes``. The generated ``Base64Str`` would reject them, since a
  CMS signature is not UTF-8 text.
* Python code works with plain values, the Base64 layer exists only on the wire (see below).
* :meth:`SignedRequest.signing_payload` and the ``verify_signature`` methods of the responses
  implement the signature rules that the test server was observed to apply.

Wire format versus Python values
--------------------------------

Every request field and the response fields ``psn`` and ``signature`` are ``format: byte``, Base64
strings on the wire. In Python they hold the plain value: text fields are ``str``, opaque fields
are ``bytes``. The Base64 layer belongs to the JSON representation only:

* ``Model.model_validate_json(body)`` reads a wire body (decodes Base64).
* ``model.model_dump_json()`` writes a wire body (encodes Base64, uses the wire field names).
* ``Model(an=..., kvnr=..., ...)`` and ``model.model_dump()`` work with plain values; the dump
  uses the wire names too.

So read an HTTP response with ``model_validate_json(response.content)``, not with
``model_validate(response.json())``: a Python ``dict`` is taken to hold plain values.

The response field ``vn`` is the exception: the documents type it as a plain string, and the VST
echoes the Vorgangsnummer exactly as the LE sent it, which is the Base64 of its hexadecimal text.
The models keep that string unchanged, because the response signature is made over it, and
offer the decoded value as :attr:`VnResult.vorgangsnummer` and :attr:`PseudonymResult.vorgangsnummer`.

Signatures
----------

A sender signs the request fields with a detached CMS signature over the plain field values
joined by ``_`` (0x5F); :meth:`SignedRequest.signing_payload` builds those bytes. The VST answers
with a SHA-256 digest over the response fields joined the same way, hashing ``psn`` decoded and
``vn`` as its wire text; ``verify_signature`` checks it. Making the CMS signature, the
Authorization token and the transport are not part of the data model.
"""

from __future__ import annotations

import base64
import binascii
import hashlib
from typing import Annotated, Any

from pydantic import (
    BaseModel,
    BeforeValidator,
    ConfigDict,
    Field,
    PlainSerializer,
    StringConstraints,
    TypeAdapter,
    ValidationInfo,
)

from ..submission.metadata import ClinicalDataNodeId, GenomicDataCenterId, SubmitterId, Tan

__all__ = [
    "Arbeitsnummer",
    "DataNodeTid",
    "ErrorDetails",
    "Kvnr",
    "MigrationRequestData",
    "Pseudonym",
    "PseudonymRequestData",
    "PseudonymResult",
    "Signature",
    "SignedRequest",
    "SubmitterTid",
    "VnRequestData",
    "VnResult",
    "Vorgangsnummer",
    "VstBaseModel",
]


# --------------------------------------------------------------------------- #
# Base64 wire layer
# --------------------------------------------------------------------------- #


def _decode_base64(value: str) -> bytes:
    try:
        return base64.b64decode(value, validate=True)
    except (binascii.Error, ValueError) as e:
        raise ValueError(f"value is not valid Base64: {e}") from e


def _decode_text(value: Any, info: ValidationInfo) -> Any:
    """In JSON mode, turn the Base64 wire value into the UTF-8 text it carries."""
    if info.mode != "json" or not isinstance(value, str):
        return value
    try:
        return _decode_base64(value).decode("utf-8")
    except UnicodeDecodeError as e:
        raise ValueError("decoded value is not UTF-8 text") from e


def _decode_bytes(value: Any, info: ValidationInfo) -> Any:
    """In JSON mode, turn the Base64 wire value into the bytes it carries."""
    if info.mode != "json" or not isinstance(value, str):
        return value
    return _decode_base64(value)


def _encode_bytes(value: bytes) -> str:
    return base64.b64encode(value).decode("ascii")


def _encode_text(value: str) -> str:
    return _encode_bytes(value.encode("utf-8"))


_DECODE_TEXT = BeforeValidator(_decode_text)
_DECODE_BYTES = BeforeValidator(_decode_bytes)
_ENCODE_TEXT = PlainSerializer(_encode_text, return_type=str, when_used="json")
_ENCODE_BYTES = PlainSerializer(_encode_bytes, return_type=str, when_used="json")


# --------------------------------------------------------------------------- #
# Field types
# --------------------------------------------------------------------------- #

Arbeitsnummer = Annotated[str, StringConstraints(min_length=3, max_length=40), _DECODE_TEXT, _ENCODE_TEXT]
"""Arbeitsnummer (AN) of the Leistungserbringer, 3 to 40 characters of text.
It must never contain patient-identifying data."""

Kvnr = Annotated[str, StringConstraints(pattern=r"^[A-Z][0-9]{9}$"), _DECODE_TEXT, _ENCODE_TEXT]
"""Unchangeable part of the Krankenversichertennummer (KVNR) per § 290 SGB V: one letter and
nine digits. The test server only knows the range ``A1111xxxxP``."""

Vorgangsnummer = Annotated[Tan, _DECODE_TEXT, _ENCODE_TEXT]
"""Vorgangsnummer (VN) of a request: 32 random bytes in their 64-character hexadecimal form, the
same value as the ``tanG``/``tanC`` of a submission. On the wire it is the Base64 of that text."""

SubmitterTid = Annotated[SubmitterId, _DECODE_TEXT, _ENCODE_TEXT]
"""Teilnehmer-ID (TID) of the Leistungserbringer: nine digits, the ``submitterId`` of a submission."""

DataNodeTid = Annotated[GenomicDataCenterId | ClinicalDataNodeId, _DECODE_TEXT, _ENCODE_TEXT]
"""Teilnehmer-ID (TID) of the data node: a GRZ-TID or a KDK-TID."""

Pseudonym = Annotated[bytes, _DECODE_BYTES, _ENCODE_BYTES]
"""Pseudonym (PSN) the VST issues for a Vorgangsnummer. Opaque; its Base64 wire form
(``model_dump(mode="json")["psn"]``) is the value to store or display."""

Signature = Annotated[bytes, _DECODE_BYTES, _ENCODE_BYTES]
"""Signature of a message: a detached CMS container from a sender, a SHA-256 digest from the VST."""

_TAN = TypeAdapter(Tan)


def _vorgangsnummer_of_response(wire: str) -> Tan:
    """Decode the ``vn`` string of a response: Base64 of the hexadecimal text as observed, or that text itself."""
    try:
        return _TAN.validate_python(wire)
    except ValueError:
        pass
    try:
        return _TAN.validate_python(_decode_base64(wire).decode("utf-8"))
    except (ValueError, UnicodeDecodeError) as e:
        raise ValueError(f"vn {wire!r} is neither a hexadecimal Vorgangsnummer nor the Base64 of one") from e


# --------------------------------------------------------------------------- #
# Models
# --------------------------------------------------------------------------- #


class VstBaseModel(BaseModel):
    """Base of all VST message bodies: wire field names as aliases, unknown fields forbidden."""

    model_config = ConfigDict(
        extra="forbid",
        validate_assignment=True,
        validate_by_name=True,
        validate_by_alias=True,
        serialize_by_alias=True,
    )


class SignedRequest(VstBaseModel):
    """A request body whose ``signature`` field covers the other fields, which subclasses declare in wire order."""

    def signing_payload(self) -> bytes:
        """
        Bytes over which the ``Signature`` of this request is made.

        It is the plain value of every field except ``Signature``, in wire order, joined by ``_``
        (0x5F): text fields as UTF-8, the Vorgangsnummer as its hexadecimal text. The order of
        the migrate and result requests was accepted by the test server on 2026-11-02; the vn
        request follows the RKI C# example.

        :returns: the payload to sign.
        """
        parts = []
        for name in type(self).model_fields:
            if name == "signature":
                continue
            value = getattr(self, name)
            parts.append(value if isinstance(value, bytes) else str(value).encode("utf-8"))
        return b"\x5f".join(parts)


class ErrorDetails(BaseModel):
    """
    Definiert einen Datensatz für die Beschreibung eines Problems

    The error body of every endpoint; the test server answers for example ``{"detail": "NotFound"}``
    for an unknown or already consumed Vorgangsnummer.
    """

    model_config = ConfigDict(
        extra="forbid",
    )
    detail: str = Field(..., min_length=1)
    """
    Eine für dieses Auftreten des Problems spezifische, für Menschen verständliche Erklärung.
    """


class VnRequestData(SignedRequest):
    """
    Definiert einen Datensatz

    Request body of ``psngen/api/v1/pseudonym/vn`` (Phase 1): the Leistungserbringer reports its
    data without a Vorgangsnummer and gets one back in a :class:`VnResult`.
    """

    an: Arbeitsnummer = Field(alias="AN")
    """
    Definiert die Arbeitsnummer des LEs.
    """
    kvnr: Kvnr = Field(alias="KVNR")
    """
    Definiert die Identifikationsnummer des Versicherten.
    """
    data_node_id: DataNodeTid = Field(alias="DataNodeId")
    """
    Definiert die Teilnehmer-Id des Knotens(GRZ-TID oder KDK-TID).
    """
    submitter_id: SubmitterTid = Field(alias="SubmitterId")
    """
    Definiert die Leistungserbringer-Teilnehmer-Id
    """
    signature: Signature = Field(alias="Signature")
    """
    Signatur des SHA256-Hash des Datenpaketes
    """


class VnResult(VstBaseModel):
    """
    Klasse, die die Vorgangsnummer kapselt

    Response body of ``psngen/api/v1/pseudonym/vn`` (Phase 1). Not yet observed from the test
    server; the field types follow the OpenAPI document and the ``vn`` handling mirrors
    :class:`PseudonymResult`.
    """

    vn: str | None = None
    """
    Gibt die Vorgangsnummer zurück, über die das Pseudonym abgerufen werden kann.
    """
    signature: Signature | None = None
    """
    Gibt die Signatur des Datenpakets zurück
    """

    @property
    def vorgangsnummer(self) -> Tan | None:
        """The issued Vorgangsnummer in its hexadecimal form, decoded from :attr:`vn`."""
        return None if self.vn is None else _vorgangsnummer_of_response(self.vn)

    def verify_signature(self) -> bool:
        """
        Check the VST's signature: the SHA-256 digest of ``vn`` as it appears in the JSON.

        By analogy with :meth:`PseudonymResult.verify_signature`; not yet confirmed against the
        test server.

        :returns: whether the signature matches the fields.
        """
        if self.signature is None or self.vn is None:
            return False
        return hashlib.sha256(self.vn.encode("utf-8")).digest() == self.signature


class MigrationRequestData(SignedRequest):
    """
    Definiert einen Datensatz

    Request body of ``psngen/api/v1/pseudonym/migrate`` (Phase 0): the Leistungserbringer reports
    its data together with the Vorgangsnummer it made itself. The test server accepts the same
    Vorgangsnummer again, which makes an already fetched Pseudonym retrievable once more.
    """

    an: Arbeitsnummer = Field(alias="AN")
    """
    Definiert die Arbeitsnummer des LEs.
    """
    kvnr: Kvnr = Field(alias="KVNR")
    """
    Definiert die Identifikationsnummer des Versicherten.
    """
    vn: Vorgangsnummer = Field(alias="VN")
    """
    Definiert die Transaktionsnummer VNg/VNk.
    """
    data_node_id: DataNodeTid = Field(alias="DataNodeId")
    """
    Definiert die Teilnehmer-Id des Knotens(GRZ-TID oder KDK-TID).
    """
    submitter_id: SubmitterTid = Field(alias="SubmitterId")
    """
    Definiert die Leistungserbringer-Teilnehmer-Id
    """
    signature: Signature = Field(alias="Signature")
    """
    Signatur des SHA256-Hash des Datenpaketes
    """


class PseudonymRequestData(SignedRequest):
    """
    Definiert einen Datensatz

    Request body of ``psn/api/v1/pseudonym/result``: a data node asks for the Pseudonym of a
    Vorgangsnummer. ``submitter_id`` is the TID of the Leistungserbringer that reported (Phase 0)
    or obtained (Phase 1) the Vorgangsnummer. A Vorgangsnummer is consumed by its first successful
    retrieval; the next request answers 404.
    """

    vn: Vorgangsnummer = Field(alias="VN")
    """
    Definiert die Transaktionsnummer VNg/VNk.
    """
    submitter_id: SubmitterTid = Field(alias="SubmitterId")
    """
    Definiert die ID des Leistungserbringer.
    """
    signature: Signature = Field(alias="Signature")
    """
    Signatur des SHA256-Hash des Datenpaketes
    """


class PseudonymResult(VstBaseModel):
    """
    Klasse, die das Pseudonym kapselt

    Response body of ``psn/api/v1/pseudonym/result``. The fields are optional in the OpenAPI
    document; the test server was observed to send all three.
    """

    psn: Pseudonym | None = None
    """
    Gibt das Pseudonym zurück.
    """
    vn: str | None = None
    """
    Gibt die Vorgangsnummer zurück, über die das Pseudonym abgerufen wurde.
    """
    signature: Signature | None = None
    """
    Gibt die Signatur des Datenpakets zurück
    """

    @property
    def vorgangsnummer(self) -> Tan | None:
        """The Vorgangsnummer in its hexadecimal form, decoded from :attr:`vn`."""
        return None if self.vn is None else _vorgangsnummer_of_response(self.vn)

    def verify_signature(self) -> bool:
        """
        Check the VST's signature.

        It is the SHA-256 digest of the response fields joined by ``_`` (0x5F): ``psn`` decoded
        (it is ``format: byte``) and ``vn`` as it appears in the JSON (it is a plain string).
        Absent fields are left out. Confirmed against the test server on 2026-11-02.

        :returns: whether the signature matches the fields.
        """
        if self.signature is None:
            return False
        parts = []
        if self.psn is not None:
            parts.append(self.psn)
        if self.vn is not None:
            parts.append(self.vn.encode("utf-8"))
        return hashlib.sha256(b"\x5f".join(parts)).digest() == self.signature
