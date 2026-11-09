# RKI Vertrauensstelle (VST) models, explained

The **Vertrauensstelle** of the Modellvorhaben Genomsequenzierung at the RKI issues the **Vorgangsnummer** (the `tanG`/`tanC` a submission carries) and the **Pseudonym** (`psn`) of a case.
Leistungserbringer (LE) report cases to it, and data nodes (GRZ, KDK) fetch the pseudonym for a Vorgangsnummer.
`grz_pydantic_models.vst` holds the request and response bodies of that REST interface.

## Where the models come from

The RKI publishes the interface as three OpenAPI documents, one per endpoint.
They are kept verbatim in `grz_pydantic_models/resources/vst/` and are the authoritative source.
The PDF specification v1.0 and the test handout describe the same interface in prose and lag behind the documents in places, for example in the casing of the `Authorization` scheme.

`vst/v1.py` started as the output of [datamodel-code-generator](https://koxudaxi.github.io/datamodel-code-generator/) over the three documents merged into one, with the command recorded in the module docstring.
It was then edited by hand.
The classes, the wire field names, which fields are required and the rejection of unknown fields are the generated ones, and the German field docstrings are the OpenAPI descriptions.
The hand edits add the field validation, the `bytes` type of `psn` and `Signature`, the plain-value Python API and the signature helpers described below.
To pick up a revised document, regenerate into a scratch file and diff it against the module.

## The three endpoints

| Endpoint (`POST`)                                                | Phase   | Sender    | Request model          | Response |
| ---------------------------------------------------------------- | ------- | --------- | ---------------------- | -------- |
| `psngen/api/v1/pseudonym/migrate`                                | 0       | LE        | `MigrationRequestData` | HTTP 200, empty body |
| `psngen/api/v1/pseudonym/vn` (`.../generate` on the test server) | 1       | LE        | `VnRequestData`        | `VnResult` |
| `psn/api/v1/pseudonym/result`                                    | 0 and 1 | GRZ / KDK | `PseudonymRequestData` | `PseudonymResult` |

In Phase 0 the LE makes the Vorgangsnummer itself and reports it afterwards (`migrate`).
In Phase 1 the LE reports the case first and the VST issues the Vorgangsnummer (`vn`).
Either way a data node later asks for the pseudonym with `result`.
Errors (400, 401, 403, 404, 415, 500) come back as an `ErrorDetails` body with a single `detail` text.

## Wire format and Python values

The VST transports the request fields and the `psn` and `signature` of a response as Base64 (OpenAPI `format: byte`).
The models keep that layer in the JSON representation only:

```python
from grz_pydantic_models.vst import PseudonymRequestData, PseudonymResult

request = PseudonymRequestData(vn=tan_g, submitter_id="260930608", signature=cms_bytes)
body = request.model_dump_json()  # Base64 values, VST field names

result = PseudonymResult.model_validate_json(response.content)
assert result.verify_signature()
psn_bytes = result.psn
psn_text = result.model_dump(mode="json")["psn"]  # the Base64 text, the form to store or show
tan_g = result.vorgangsnummer  # the 64 hexadecimal characters
```

Read responses with `model_validate_json`, not `model_validate(response.json())`: a Python `dict` is taken to hold plain values.
`model_dump()` gives the plain values under the wire names.

| Field          | Wire name                | Python                    | Reuses |
| -------------- | ------------------------ | ------------------------- | ------ |
| Arbeitsnummer  | `AN`                     | `str`, 3 to 40 characters |        |
| KVNR           | `KVNR`                   | `str`, `^[A-Z][0-9]{9}$`  |        |
| Vorgangsnummer | `VN`                     | `str`, 64 hex characters  | `Tan` of the submission metadata |
| LE-TID         | `SubmitterId`            | `str`, 9 digits           | `SubmitterId` of the submission metadata |
| GRZ-/KDK-TID   | `DataNodeId`             | `str`                     | `GenomicDataCenterId \| ClinicalDataNodeId` |
| Pseudonym      | `psn`                    | `bytes`                   |        |
| Signature      | `Signature`, `signature` | `bytes`                   |        |

Requests use PascalCase names, responses lowercase ones, exactly as the documents say.
Unknown fields raise, as they do at the VST.

## The Vorgangsnummer on the wire

A request carries the Vorgangsnummer as the Base64 of its 64-character hexadecimal text, which is what the handout defines and the RKI C# examples send.
The PDF's worked examples carry the Base64 of the 32 raw bytes instead; the models reject that form.

A response types `vn` as a plain string, not as `format: byte`.
The test server echoes the string the LE sent, so in practice it is the same Base64 text.
The models keep `vn` exactly as received, because the VST signs that text, and decode it on demand as `vorgangsnummer`, which also accepts the bare hexadecimal text.

## Signatures

`signing_payload()` of a request gives the bytes the sender signs with a detached CMS/PKCS#7 signature: the plain value of every field but `Signature`, in wire order, joined by `_` (0x5F).
The migrate and result orders were accepted by the test server on 2026-11-02; the vn request follows the RKI C# example.

`verify_signature()` of a response checks the VST's signature, a SHA-256 digest over the response fields joined by `_`: `psn` decoded, since it is `format: byte`, and `vn` as it appears in the JSON, since it is a plain string.
The formula was confirmed against three server-signed responses of 2026-11-02, which are the fixtures in `tests/data/vst/`.
For `VnResult` the analogous digest over `vn` alone is implemented but unconfirmed, since Phase 1 has not been exercised yet.

Making the CMS signature, the `Authorization: Custom <token>` header and the HTTP transport belong to a client, not to the data model.

## Observed behaviour of the test server

These points are not in the RKI documents and matter for a client:

- A Vorgangsnummer is consumed by its first successful `result` call; the same request answers `404 {"detail": "NotFound"}` afterwards.
- `migrate` accepts the same Vorgangsnummer again, after which `result` returns the same pseudonym once more.
  So after a lost `result` response the GRZ asks the LE to re-deliver the Vorgangsnummer; it cannot recover alone.
- The pseudonym is stable per KVNR: different Vorgangsnummern for the same KVNR returned the same `psn`.
- The `Authorization` scheme is case-sensitive: `Custom` passes, `custom` answers a bare `401 {"detail": "Unauthorized"}`.
- The Phase 1 endpoint is served as `.../pseudonym/generate`, while the OpenAPI document names it `.../pseudonym/vn`.

## Open points against the live VST

- Whether the response `vn` is always the LE's string echoed verbatim, and what `VnResult` returns for a Vorgangsnummer the VST issued itself.
- The `VnResult` signature formula.
- The OpenAPI document of the `vn` endpoint declares its 200 response without a body, although it defines `VnResult`.
