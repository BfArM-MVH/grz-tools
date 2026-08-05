# MII consent terminology and profile, vendored verbatim

Conformance inputs for the `grz_pydantic_models.mii.consent` model, copied unchanged from the
published FHIR packages of `de.medizininformatikinitiative.kerndatensatz.consent`.

The consent tests **discover** these files. Adding one adds test cases; nothing else has to be
registered anywhere.

## Adding a newly released consent package

```bash
python packages/grz-pydantic-models-testing/scripts/vendor_mii_consent_package.py 2026.0.1
uv run pytest packages/grz-pydantic-models
```

The script fetches the package, writes each artefact under the name given below, and reports whether
it was added, was already vendored unchanged, or is not shipped by that package. It refuses to
overwrite a vendored file whose contents differ from upstream and exits non-zero instead: the MII has
published two different profiles under version 1.0.9, so a version number alone is not proof that
nothing changed. Pass `--force` once you have looked at the difference and want upstream to win.

Only a human runs this. The test suite never reaches the network, so a registry outage cannot turn a
build red and every change to these artefacts arrives as a reviewable diff.

### Doing it by hand

```bash
VERSION=2026.0.1   # the package version the MII published
curl -sSL -o pkg.tgz "https://packages2.fhir.org/packages/de.medizininformatikinitiative.kerndatensatz.consent/$VERSION"
tar xzf pkg.tgz
```

From `package/`, copy each of the three artefacts below into this directory, naming it after the
**version stated inside the resource**, not after the package version:

| Artefact in `package/`                              | File name here                                                   |
| --------------------------------------------------- | ---------------------------------------------------------------- |
| `CodeSystem-MiiConsentPolicyCodeSystem.json`        | `mii-cs-consent-policy.<CodeSystem.version>.json`                |
| `CodeSystem-MiiConsentVersionModuleCodeSystem.json` | `mii-cs-consent-version-modules.<CodeSystem.version>.json`       |
| `Profile_MII_Consent_Einwilligung.json`             | `mii-pr-consent-einwilligung.<StructureDefinition.version>.json` |

Several packages ship the same resource version, so a file may already be present and identical; in
that case there is nothing to add. Copy the JSON verbatim apart from re-indenting.

Then run the consent tests. Anything the model does not yet cover is reported by file name:

- an OID added to the version and module CodeSystem that `BROAD_CONSENT_DOCUMENT_OIDS` does not
  classify
- either policy code research consent is derived from going missing, deprecated or inactive
- a profile that pins a different category CodeSystem, relaxes or tightens a period bound, changes
  the minimum number of categories, or stops forbidding what the model rejects

A failure is a prompt to update the model, not a reason to edit these files.

## Which package ships which version

| Package  | Policy CodeSystem | Version and modules CodeSystem | Profile |
| -------- | ----------------- | ------------------------------ | ------- |
| 2025.0.1 | 1.0.5             | not shipped                    | 1.0.8   |
| 2025.0.2 | 1.0.6             | not shipped                    | 1.0.8   |
| 2025.0.3 | 1.0.7             | not shipped                    | 1.0.8   |
| 2025.0.4 | 1.0.7             | not shipped                    | 1.0.8   |
| 2026.0.0 | 1.1.0             | 0.2.0                          | 1.0.9   |

Package 2026.0.0 is the only published package defining the version and module CodeSystem. The 2025
packages use its OIDs as `Consent.policy[].uri` without defining them, and the 2026.0.1 release
candidates drop the file again while still pinning its canonical URL in the profile.
