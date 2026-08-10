# MII research consent, explained

A donor signs a consent form on paper: the **MII broad consent** ("breite Einwilligung" of the
Medizininformatik-Initiative). The hospital records that decision in its consent management
system, which exports it as a **FHIR `Consent` resource**. A GRZ submission carries that resource
inside `researchConsents[].scope`, and this package answers two questions about it:

1. **Is it well-formed?** Does it follow the MII's FHIR specification (validation)?
2. **What does it say?** May the donor's data be used for research at a given date (evaluation)?

Current version numbers live in the code constants and in the `example_terminology` README of
`grz-pydantic-models-testing`.

## Glossary

| Term | In one sentence |
| --- | --- |
| MII broad consent | The standardized consent form a donor signs, allowing broad research use of their data. Exists in several document versions. |
| KDS consent package | The MII's FHIR package (`de.medizininformatikinitiative.kerndatensatz.consent`) that defines how to represent such a consent digitally ([registry](https://packages2.fhir.org/packages/de.medizininformatikinitiative.kerndatensatz.consent)). |
| Profile | The rulebook inside the package (`MII_PR_Consent_Einwilligung`): which fields a Consent resource must/may have. |
| Policy CodeSystem | The catalogue of **permissions** ("use data scientifically", "contact again", …), one OID each, organized as modules with sub-items. |
| Version & modules CodeSystem | The catalogue of **signable documents** (each broad consent version, withdrawal forms, additional modules), one OID each. |
| GRZ metadata schema | BfArM's format for the whole submission; `researchConsents[]` is one small part of it. |

## A consent, annotated

```jsonc
{
  "resourceType": "Consent",
  "status": "active",              // modifier: anything else means "not in force"
  "scope": { "coding": [ { "system": ".../consentscope", "code": "research" } ] },
  "category": [
    { "coding": [ { "system": "http://loinc.org", "code": "57016-8" } ] },
    { "coding": [ { /* the MII broad consent category */ } ] }
  ],
  "patient": { "reference": "Patient/..." },     // or an identifier with system + value
  "dateTime": "2020-09-01",
  "policy": [
    { "uri": "urn:oid:2.16.840.1.113883.3.1937.777.24.2.184" }
    // ^ WHICH document was signed (an OID from the version & modules CodeSystem)
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

- `policy[].uri` says **which document** the donor signed → version & modules CodeSystem.
- `provision.provision[].code` says **what the donor allowed or refused** → policy CodeSystem.

## Four version numbers that mean different things

| Version of… | Where you see it | Source of truth in this repo |
| --- | --- | --- |
| the GRZ metadata schema | `$schema` URL of the submission | `is_supported_version`, `get_accepted_versions` |
| the KDS consent package | `researchConsents[].schemaVersion` | `RESEARCH_CONSENT_SCHEMA_VERSIONS` |
| each artefact inside the package (profile, CodeSystems) | inside the package files | `example_terminology/` file names + its README's mapping table |
| the signed broad consent document | OIDs in `Consent.policy[].uri` | `BroadConsentVersion`, `BROAD_CONSENT_DOCUMENT_OIDS` |

None of these implies another. The classic trap: **`schemaVersion` names the KDS package (the
digital format), not the document the donor signed.** Which broad consent version the donor
actually signed is derived from the policy OIDs (`Consent.broad_consent_versions`).

Two facts that follow from the axes moving independently:

- Most package releases only grow the terminology and leave the profile untouched, so one and the
  same Consent stays valid under every accepted `schemaVersion` (a test asserts this).
- The GRZ metadata schema sometimes lists a package version before the MII has released it; such
  versions stay rejected until the release exists and its artefacts are vendored.

## Question 1: is it well-formed?

The model (`grz_pydantic_models/mii/consent.py`) enforces what the profile pins:

| Element | Rule |
| --- | --- |
| `status` | required; only `active` marks the consent as in force |
| `scope` | exactly one coding, and it must be `research` |
| `category` | must contain the LOINC consent category and the MII broad consent category, one coding each; extra categories are allowed (open slicing) |
| `patient` | must identify the patient: reference or identifier (identifier needs `system` + `value`) |
| `policy` | at least one document OID, with or without the `urn:oid:` prefix |
| `provision` (root) | `type` must be `deny` (opt-in), `period` required, `code` forbidden |
| `provision.provision[]` | the decisions: `type`, `period`, at least one `code` |
| a third provision level | forbidden |

Two quirks worth knowing:

- The MII renamed the category CodeSystem at one point but kept shipping examples in the old
  spelling, so **both spellings are in the wild and both are accepted**
  (`MII_CONSENT_CATEGORY_SYSTEMS`).
- Profile 1.0.9, shipped by package 2026.0.0, allows a `period` without an `end`, and such a
  period never expires. Profile 1.0.8, shipped by every 2025 package, still requires an `end`, so
  an open-ended period is **rejected under those `schemaVersion`s** and accepted only under
  2026.0.0 (`PROFILES_REQUIRING_PERIOD_END`). The field is optional on `Period` itself, because
  one model serves every profile version; the version that decides is the one the submission
  declares.

Unknown fields are generally ignored (FHIR resources carry much more than we read), **except**
where ignoring would silently drop a permission the evaluation never looks at: a `code` on the
root provision, a third provision level, a wrong `resourceType`, an unidentified patient. Those
are rejected so the submitter notices.

## Question 1b: what was signed?

Every OID the MII ships in the version & modules CodeSystem is classified in
`BROAD_CONSENT_DOCUMENT_OIDS` as one of:

- **consent**: the broad consent itself, per document version (incl. minor / legal-guardian variants)
- **rejection**: the donor refused (Ablehnung)
- **complete withdrawal**: everything revoked (Komplettwiderruf)
- **partial withdrawal**: parts revoked (Teilwiderruf)
- **additional module**: an add-on to the broad consent (Zusatzmodul)

An OID we do not recognize (e.g. a future broad consent version) is **reported, never rejected**
(`unknown_document_oids`): new document versions must not break submissions.

## Question 2: does the donor consent to research?

Research use requires one specific permission from the policy CodeSystem, and that catalogue is
hierarchical: permissions live inside modules:

```text
…24.5.3.1  "PATDAT erheben, speichern, nutzen"      ← the module
└── …24.5.3.8  "MDAT wissenschaftlich nutzen"       ← the permission research needs
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
in `example_terminology/` of `grz-pydantic-models-testing`. Its README describes the
package-to-artefact mapping and the update workflow.
