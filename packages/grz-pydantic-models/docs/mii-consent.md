# MII research consent, explained

A donor signs a consent form on paper: the **MII broad consent** ("breite Einwilligung" of the
Medizininformatik-Initiative). The hospital records that decision in its consent management
system, which exports it as a **FHIR `Consent` resource**. A GRZ submission carries that resource
inside `researchConsents[].scope`, and this package answers two questions about it:

1. **Is it well-formed?** Does it follow the MII's FHIR specification (validation)?
2. **What does it say?** May the donor's data be used for research at a given date (evaluation)?

`scope` may also be absent, paired with a `noScopeJustification` saying why. And a `scope` that
fails Consent validation is not rejected: it is kept as a plain object, the validation errors are
logged as warnings, and both questions below then answer "nothing known" rather than "no".

Current version numbers live in the code constants and in `example_terminology/packages.json` of
`grz-pydantic-models-testing`.

## Glossary

| Term | In one sentence |
| --- | --- |
| MII broad consent | The standardized consent form a donor signs, allowing broad research use of their data. Exists in several document versions. |
| KDS consent package | The MII's FHIR package (`de.medizininformatikinitiative.kerndatensatz.consent`) that defines how to represent such a consent digitally ([registry](https://packages2.fhir.org/packages/de.medizininformatikinitiative.kerndatensatz.consent)). |
| Profile | The rulebook inside the package (`MII_PR_Consent_Einwilligung`): which fields a Consent resource must/may have. |
| Policy CodeSystem | The catalogue of **permissions** ("use data scientifically", "contact again", …), one OID each, organized as modules with sub-items. |
| Version and modules CodeSystem | The catalogue of **signable documents** (each broad consent version, withdrawal forms, additional modules), one OID each. |
| GRZ metadata schema | BfArM's format for the whole submission; `researchConsents[]` is one small part of it. |

## A consent, annotated

```jsonc
{
  "resourceType": "Consent",
  "status": "active",              // modifier: anything else means "not in force"
  "scope": { "coding": [ { "system": ".../consentscope", "code": "research" } ] },
  "category": [
    { "coding": [ { "system": "http://loinc.org", "code": "57016-8" } ] },
    { "coding": [ { "system": ".../mii-cs-consent-consent_category",
                    "code": "2.16.840.1.113883.3.1937.777.24.2.184" } ] }
    // ^ THAT this is an MII broad consent at all, not which version of it
  ],
  "patient": { "reference": "Patient/..." },     // or an identifier with system + value
  "dateTime": "2020-09-01",
  "policy": [
    { "uri": "urn:oid:2.16.840.1.113883.3.1937.777.24.2.1791" }
    // ^ WHICH document was signed: broad consent 1.6f, an OID from the version and modules CodeSystem
  ],
  "provision": {                   // root provision: the opt-in frame
    "type": "deny",                // nothing is allowed unless a nested rule permits it
    "period": { "start": "2020-09-01", "end": "2050-08-31" },   // bounds ALL nested rules
    "provision": [                 // the actual decisions, one per permission (group)
      {
        "type": "permit",
        "period": { "start": "2020-09-01", "end": "2025-08-31" },
        "code": [ { "coding": [ {
          "code": "2.16.840.1.113883.3.1937.777.24.5.3.1"
          // ^ WHAT is permitted (an OID from the policy CodeSystem)
        } ] } ]
      }
    ]
  }
}
```

Two different OID catalogues meet here, and telling them apart is the key to everything else:

- `policy[].uri` says **which document** the donor signed → version and modules CodeSystem.
- `provision.provision[].code` says **what the donor allowed or refused** → policy CodeSystem.

Some systems state the signed document as a `category` coding from the version and modules
CodeSystem instead of in `policy[]`. `Consent.document_oids` reads both places, so either spelling
works.

## Four version numbers that mean different things

| Version of… | Where you see it | Source of truth in this repo |
| --- | --- | --- |
| the GRZ metadata schema | `$schema` URL of the submission | `is_supported_version`, `get_accepted_versions` |
| the KDS consent package | `researchConsents[].schemaVersion` | `RESEARCH_CONSENT_SCHEMA_VERSIONS` |
| each artefact inside the package (profile, CodeSystems) | inside the package files | `example_terminology/` file names + `example_terminology/packages.json` |
| the signed broad consent document | OIDs in `Consent.policy[].uri`, or in a `category` coding from the version and modules CodeSystem | `BroadConsentVersion`, `BROAD_CONSENT_DOCUMENT_OIDS` |

None of these implies another. The classic trap: **`schemaVersion` names the KDS package (the
digital format), not the document the donor signed.** Which broad consent version the donor
actually signed is derived from those document OIDs (`Consent.broad_consent_versions`), wherever
they appear.

Two facts that follow from the axes moving independently:

- Most package releases only grow the terminology and leave the profile untouched, so a Consent
  that parses under one `schemaVersion` usually parses under all of them
  (`test_every_published_package_version_parses_the_same_consent` pins this for a fully bounded
  consent). Where the profile did change, portability ends: see the period rule below.
- The GRZ metadata schema sometimes lists a package version before the MII has released it; such
  versions stay rejected until the release exists and its artefacts are vendored.

## Question 1: is it well-formed?

The model (`grz_pydantic_models/mii/consent.py`) enforces what the profile pins, plus two rules the
profile leaves open but the GRZ requires, marked †. The root-provision rule is enforced one level
up, on `ResearchConsent` in `submission/metadata/v1.py`.

| Element | Rule |
| --- | --- |
| `status` | required; only `active` marks the consent as in force |
| `scope` | exactly one coding, and it must be `research` |
| `category` | must contain the LOINC consent category and the MII broad consent category, one coding each; extra categories are allowed (open slicing) |
| `patient` | † must identify the patient: reference or identifier (identifier needs `system` + `value`). The profile requires `patient` but marks both ways of filling it mustSupport only |
| `policy` | at least one document OID, with or without the `urn:oid:` prefix |
| `provision` (root) | † `type` must be `deny` (opt-in); the profile requires `type` but fixes no value. `period` required, `code` forbidden |
| `provision.provision[]` | the decisions: `type`, `period`, at least one `code` |
| a third provision level | forbidden |

Two quirks worth knowing:

- The MII renamed the category CodeSystem at one point but kept shipping examples in the old
  spelling, so **both spellings are in the wild and both are accepted**
  (`MII_CONSENT_CATEGORY_SYSTEMS`).
- Profile 1.0.9, shipped by package 2026.0.0, relaxed both provision `period.end` elements from
  1..1 to 0..1. FHIR reads a period with no `end` as still running, so such a permission never
  expires. Profile 1.0.8, shipped by every 2025 package, still requires an `end`, so an open-ended
  period is **rejected under those `schemaVersion`s** (`PROFILES_REQUIRING_PERIOD_END`). It is
  accepted under 2026.0.0, and also when the submission declares no `schemaVersion` at all, which
  metadata 1.3 permits: with no package named there is no profile to enforce. The field stays
  optional on `Period` itself, because one model serves every profile version; the version that
  decides is the one the submission declares.

Unknown fields are ignored on the resource and its provisions (FHIR carries much more than we
read), but **forbidden** on the small leaf elements the evaluation reads field by field: `Period`,
`Coding` and `CodeableConcept`. Four things are rejected outright rather than ignored: a `code` on
the root provision and a third provision level, because they carry permissions the evaluation never
looks at and accepting them would silently lose them; a wrong `resourceType` and a patient
identified by neither reference nor identifier, because neither can be what the submitter meant.

## Question 1b: what was signed?

Every OID the MII ships in the version and modules CodeSystem is classified in
`BROAD_CONSENT_DOCUMENT_OIDS` as one of:

- **consent**: the broad consent itself, per document version (incl. minor / legal-guardian variants)
- **rejection**: the donor refused (Ablehnung)
- **complete withdrawal**: everything revoked (Komplettwiderruf)
- **partial withdrawal**: parts revoked (Teilwiderruf)
- **additional module**: an add-on to the broad consent (Zusatzmodul)

An OID we do not recognize (e.g. a future broad consent version) is **reported, never rejected**
(`unknown_document_oids`): new document versions must not break submissions.

**None of this feeds the research decision.** Question 2 reads only `status` and the provisions; the
document kind is reported, not enforced. A real withdrawal or rejection says so twice, through its
document OID *and* by denying the permissions or stating none, and it is the second half that
Question 2 acts on. The shipped `withdrawal_complete` example shows the shape: its status is
`active`, and what refuses is its nested `deny` of the PATDAT module. So treat a withdrawal or
rejection OID on a consent that still permits research as a contradiction worth investigating, not
as a refusal the model has already applied.

## Question 2: does the donor consent to research?

Research use requires one specific permission from the policy CodeSystem, and that catalogue is
hierarchical: permissions live inside modules:

```text
…24.5.3.1  "Patientendaten erheben, speichern, nutzen"   ← the module (PATDAT_ERHEBEN_SPEICHERN_NUTZEN)
└── …24.5.3.8  "MDAT wissenschaftlich nutzen"            ← what research needs
                                                            (MDAT_WISSENSCHAFTLICH_NUTZEN_EU_DSGVO_NIVEAU;
                                                             policy 1.0.5-1.0.7 display it as "… EU DSGVO NIVEAU")
```

Permitting the module includes everything inside it, so **either** OID being permitted grants
research consent. That is why `ResearchConsentCodes` contains both.

`consents_to_research(consents, date)` then applies these rules:

1. A consent whose `status` is not `active` contributes nothing. `status` is a FHIR modifier
   element, so a rejected, withdrawn (`inactive`), draft, proposed or erroneous consent grants
   nothing regardless of its provisions.
2. A consent whose root provision period does not contain `date` contributes nothing, since the root
   deny frame bounds every nested rule. Same for each nested provision's own period.
3. Codes are matched on the OID alone, because the surrounding `system` is written
   inconsistently in the wild (with and without `urn:oid:`).
4. **Deny wins.** A permit of either research code grants, but a deny on either revokes, also
   across multiple consents of the same donor. So a partial withdrawal that denies `…5.3.8`
   revokes research use even while the module permit still stands.
5. Date-only period bounds cover the whole day (a start begins at midnight, an end expires at the
   end of that day); datetimes without a timezone are read as UTC.
6. A consent that states neither research code grants nothing: silence is not consent.

## How this stays correct over time

Every assumption above (the OID classification, the research codes' existence and hierarchy, the
profile's cardinalities and prohibitions) is re-checked against the MII's own artefacts, vendored
in `example_terminology/` of `grz-pydantic-models-testing`. Its `packages.json` records which
package ships which artefact version, and its README describes the update workflow.
