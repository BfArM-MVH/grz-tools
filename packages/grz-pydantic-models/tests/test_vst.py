"""
Tests of the RKI Vertrauensstelle (VST) models.

The fixtures in ``tests/data/vst/`` are the bodies the RKI test server exchanged with the journey
script on 2026-11-02: a migrate request, a result request and the result responses it signed.
"""

import base64
import hashlib
import json
from pathlib import Path

import pytest
from grz_pydantic_models.vst import (
    ErrorDetails,
    MigrationRequestData,
    PseudonymRequestData,
    PseudonymResult,
    VnRequestData,
    VnResult,
)
from pydantic import ValidationError

DATA = Path(__file__).parent / "data" / "vst"


def b64(raw: bytes) -> str:
    return base64.b64encode(raw).decode("ascii")


def load(name: str):
    return json.loads((DATA / name).read_text())


# Plain values of the recorded run: the fixed test values of the RKI handout and the VN of that run.
AN_TEXT = "AN-TEST-0001"
KVNR_TEXT = "A111100008"
VN_HEX = "ee9464566fd78a35eaf7908a2a27c551061942db204e2bccf5970774bef6f944"
VN_WIRE = b64(VN_HEX.encode())
DATA_NODE_ID_TEXT = "GRZK00001"
SUBMITTER_ID_TEXT = "260930608"

MIGRATE_WIRE = load("migrate_request.json")
RESULT_REQUEST_WIRE = load("result_request.json")
RESULT_RESPONSES_WIRE = load("result_responses.json")


# --------------------------------------------------------------------------- #
# Requests
# --------------------------------------------------------------------------- #


def test_recorded_migrate_request_reads_and_writes_the_wire_form():
    model = MigrationRequestData.model_validate_json(json.dumps(MIGRATE_WIRE))

    assert model.an == AN_TEXT
    assert model.kvnr == KVNR_TEXT
    assert model.vn == VN_HEX
    assert model.data_node_id == DATA_NODE_ID_TEXT
    assert model.submitter_id == SUBMITTER_ID_TEXT
    # A detached CMS container: DER SEQUENCE, not UTF-8 text.
    assert model.signature[:2] == b"\x30\x82"

    assert json.loads(model.model_dump_json()) == MIGRATE_WIRE


def test_recorded_result_request_round_trips():
    model = PseudonymRequestData.model_validate_json(json.dumps(RESULT_REQUEST_WIRE))

    assert model.vn == VN_HEX
    assert model.submitter_id == SUBMITTER_ID_TEXT
    assert json.loads(model.model_dump_json()) == RESULT_REQUEST_WIRE


def test_request_built_from_plain_values_matches_the_recorded_body():
    model = MigrationRequestData(
        an=AN_TEXT,
        kvnr=KVNR_TEXT,
        vn=VN_HEX,
        data_node_id=DATA_NODE_ID_TEXT,
        submitter_id=SUBMITTER_ID_TEXT,
        signature=base64.b64decode(MIGRATE_WIRE["Signature"]),
    )

    assert json.loads(model.model_dump_json()) == MIGRATE_WIRE
    assert model.model_dump() == {
        "AN": AN_TEXT,
        "KVNR": KVNR_TEXT,
        "VN": VN_HEX,
        "DataNodeId": DATA_NODE_ID_TEXT,
        "SubmitterId": SUBMITTER_ID_TEXT,
        "Signature": base64.b64decode(MIGRATE_WIRE["Signature"]),
    }


def test_wire_field_order_is_the_documented_one():
    assert list(json.loads(MigrationRequestData.model_validate_json(json.dumps(MIGRATE_WIRE)).model_dump_json())) == [
        "AN",
        "KVNR",
        "VN",
        "DataNodeId",
        "SubmitterId",
        "Signature",
    ]


def test_signing_payload_of_the_requests():
    """Plain values joined by ``_``: AN_KVNR_VN_DataNodeId_SubmitterId, AN_KVNR_DataNodeId_SubmitterId, VN_SubmitterId."""
    migrate = MigrationRequestData.model_validate_json(json.dumps(MIGRATE_WIRE))
    assert migrate.signing_payload() == b"_".join(
        [AN_TEXT.encode(), KVNR_TEXT.encode(), VN_HEX.encode(), DATA_NODE_ID_TEXT.encode(), SUBMITTER_ID_TEXT.encode()]
    )

    result = PseudonymRequestData.model_validate_json(json.dumps(RESULT_REQUEST_WIRE))
    assert result.signing_payload() == VN_HEX.encode() + b"_" + SUBMITTER_ID_TEXT.encode()

    vn_request = VnRequestData.model_validate_json(json.dumps({k: v for k, v in MIGRATE_WIRE.items() if k != "VN"}))
    assert vn_request.signing_payload() == b"_".join(
        [AN_TEXT.encode(), KVNR_TEXT.encode(), DATA_NODE_ID_TEXT.encode(), SUBMITTER_ID_TEXT.encode()]
    )


def test_vn_request_has_no_vn_field():
    with pytest.raises(ValidationError, match="extra_forbidden"):
        VnRequestData.model_validate_json(json.dumps(MIGRATE_WIRE))


def test_wire_names_are_case_sensitive():
    """The OpenAPI documents give exact names and the server was seen to use them; other casings are unknown fields."""
    wire = {**{k: v for k, v in MIGRATE_WIRE.items() if k != "KVNR"}, "Kvnr": MIGRATE_WIRE["KVNR"]}
    with pytest.raises(ValidationError, match="extra_forbidden"):
        MigrationRequestData.model_validate_json(json.dumps(wire))


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("AN", b64(b"ab")),  # too short
        ("AN", b64(b"x" * 41)),  # too long
        ("KVNR", b64(b"not-a-kvnr")),
        ("KVNR", b64(b"A11110000")),  # nine characters
        ("DataNodeId", b64(b" GRZK00001")),  # stray space, as in one example of the PDF spec
        ("DataNodeId", b64(b"XYZK00007")),  # neither GRZ nor KDK
        ("SubmitterId", b64(b"ID-LE-12341")),  # an example of the PDF spec, not a nine-digit TID
        ("VN", b64(b"too short")),
        ("VN", b64(b"g" * 64)),  # 64 characters, not hexadecimal
        ("VN", b64(bytes.fromhex(VN_HEX))),  # the 32 raw bytes, as in the PDF's examples, not the hex text
    ],
)
def test_invalid_field_values_are_rejected(field: str, value: str):
    wire = {**MIGRATE_WIRE, field: value}
    with pytest.raises(ValidationError):
        MigrationRequestData.model_validate_json(json.dumps(wire))


@pytest.mark.parametrize("field", ["AN", "VN", "Signature"])
def test_non_base64_value_is_rejected(field: str):
    wire = {**MIGRATE_WIRE, field: "not base64 !!!"}
    with pytest.raises(ValidationError, match="not valid Base64"):
        MigrationRequestData.model_validate_json(json.dumps(wire))


def test_non_utf8_text_is_rejected():
    wire = {**MIGRATE_WIRE, "AN": b64(b"\xff\xfe\xfd")}
    with pytest.raises(ValidationError, match="not UTF-8"):
        MigrationRequestData.model_validate_json(json.dumps(wire))


# --------------------------------------------------------------------------- #
# Responses
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("wire", RESULT_RESPONSES_WIRE, ids=lambda w: w["vn"][:12])
def test_recorded_result_responses_verify(wire: dict):
    model = PseudonymResult.model_validate_json(json.dumps(wire))

    assert len(model.psn) == 32
    assert model.vn == wire["vn"]  # kept as sent
    assert model.vorgangsnummer == base64.b64decode(wire["vn"]).decode()
    assert model.verify_signature()
    assert json.loads(model.model_dump_json()) == wire


def test_tampered_result_response_fails_verification():
    wire = {**RESULT_RESPONSES_WIRE[0], "psn": b64(bytes(32))}
    assert not PseudonymResult.model_validate_json(json.dumps(wire)).verify_signature()


def test_response_signature_hashes_the_wire_text_of_vn():
    """The hex text is what the LE signs; the VST hashes the Base64 text it echoes."""
    psn = bytes(range(32))
    over_wire_text = hashlib.sha256(psn + b"_" + VN_WIRE.encode()).digest()
    over_hex_text = hashlib.sha256(psn + b"_" + VN_HEX.encode()).digest()

    assert PseudonymResult(psn=psn, vn=VN_WIRE, signature=over_wire_text).verify_signature()
    assert not PseudonymResult(psn=psn, vn=VN_WIRE, signature=over_hex_text).verify_signature()


def test_response_fields_are_optional():
    empty = PseudonymResult.model_validate_json("{}")
    assert (empty.psn, empty.vn, empty.signature) == (None, None, None)
    assert empty.vorgangsnummer is None
    assert not empty.verify_signature()

    nulls = PseudonymResult.model_validate_json('{"psn": null, "vn": null, "signature": null}')
    assert nulls == empty


def test_response_vn_accepts_the_hexadecimal_text_too():
    assert PseudonymResult(vn=VN_HEX).vorgangsnummer == VN_HEX
    with pytest.raises(ValueError, match="neither a hexadecimal Vorgangsnummer"):
        _ = PseudonymResult(vn="bm90IGEgdm4=").vorgangsnummer  # "not a vn"


def test_vn_result_mirrors_pseudonym_result():
    signature = hashlib.sha256(VN_WIRE.encode()).digest()
    model = VnResult.model_validate_json(json.dumps({"vn": VN_WIRE, "signature": b64(signature)}))

    assert model.vorgangsnummer == VN_HEX
    assert model.verify_signature()
    assert not VnResult(vn=VN_WIRE).verify_signature()
    assert VnResult.model_validate_json("{}").vorgangsnummer is None


def test_unknown_field_is_rejected():
    with pytest.raises(ValidationError, match="extra_forbidden"):
        PseudonymResult.model_validate_json(json.dumps({**RESULT_RESPONSES_WIRE[0], "traceId": "x"}))


# --------------------------------------------------------------------------- #
# Errors
# --------------------------------------------------------------------------- #


def test_error_details():
    assert ErrorDetails.model_validate_json('{"detail": "NotFound"}').detail == "NotFound"
    for body in ("{}", '{"detail": ""}', '{"detail": "Unauthorized", "traceId": "abc"}'):
        with pytest.raises(ValidationError):
            ErrorDetails.model_validate_json(body)
