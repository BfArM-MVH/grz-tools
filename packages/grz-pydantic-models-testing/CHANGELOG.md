# Changelog

## [1.0.0](https://github.com/BfArM-MVH/grz-tools/compare/grz-pydantic-models-testing-v0.1.0...grz-pydantic-models-testing-v1.0.0) (2026-08-13)


### ⚠ BREAKING CHANGES

* **grz-pydantic-models:** a consent dateTime stating a time must spell out seconds and separate date from time by T, and a date must be a dashed calendar date, as the FHIR dateTime regex requires. Metadata stating a bound in any other form is now rejected rather than parsed. The four fields hold a FhirDateTime rather than a datetime, keeping the value as submitted and resolving to a moment through .first_moment or .last_moment. A naive datetime handed to one of them in code is refused rather than read as UTC.
* **grz-pydantic-models:** ResearchConsent.consents_to_research, and with it the consented flag grz-db records for a submission, now returns False for every consent whose status is not active. Metadata that declares a 2025 researchConsent schemaVersion and leaves any provision period open-ended is now rejected rather than accepted.

### Features

* **grz-common:** report a consent dateTime without a timezone when ([c6c4f2f](https://github.com/BfArM-MVH/grz-tools/commit/c6c4f2ff051f581430653b2efcc4b156fe5ea33f))
* **grz-db:** let a submission diff withhold overwrites outside an ([c6c4f2f](https://github.com/BfArM-MVH/grz-tools/commit/c6c4f2ff051f581430653b2efcc4b156fe5ea33f))
* **grz-pydantic-models-testing:** vendor MII consent artefacts and add ([1390ab8](https://github.com/BfArM-MVH/grz-tools/commit/1390ab8eae05452283d7a4ad6f4a1c1501bf5ed4))
* **grz-pydantic-models:** accept every released MII consent package ([#640](https://github.com/BfArM-MVH/grz-tools/issues/640)) ([1390ab8](https://github.com/BfArM-MVH/grz-tools/commit/1390ab8eae05452283d7a4ad6f4a1c1501bf5ed4)), closes [#639](https://github.com/BfArM-MVH/grz-tools/issues/639)
* **grzctl:** overwrite selected fields during db backfill ([c6c4f2f](https://github.com/BfArM-MVH/grz-tools/commit/c6c4f2ff051f581430653b2efcc4b156fe5ea33f))


### Bug Fixes

* **grz-common,grz-db:** require grz-pydantic-models 3 ([0d17b7e](https://github.com/BfArM-MVH/grz-tools/commit/0d17b7e81ded026ee41762d02b0f84c4a1b2f9b3))
* **grz-pydantic-models-testing:** require grz-pydantic-models 3 ([0d17b7e](https://github.com/BfArM-MVH/grz-tools/commit/0d17b7e81ded026ee41762d02b0f84c4a1b2f9b3))
* **grz-pydantic-models:** Change dateTime field type from date to datetime in Consent model ([#611](https://github.com/BfArM-MVH/grz-tools/issues/611)) ([2c69f9e](https://github.com/BfArM-MVH/grz-tools/commit/2c69f9e0d29e7633100b40a60d762859c9bab0a7))
* **grz-pydantic-models:** compare metadata schema versions correctly ([1390ab8](https://github.com/BfArM-MVH/grz-tools/commit/1390ab8eae05452283d7a4ad6f4a1c1501bf5ed4))
* **grz-pydantic-models:** follow the MII consent profile more closely ([1390ab8](https://github.com/BfArM-MVH/grz-tools/commit/1390ab8eae05452283d7a4ad6f4a1c1501bf5ed4))
* **grz-pydantic-models:** grant research consent only while in force ([#642](https://github.com/BfArM-MVH/grz-tools/issues/642)) ([0d17b7e](https://github.com/BfArM-MVH/grz-tools/commit/0d17b7e81ded026ee41762d02b0f84c4a1b2f9b3))
* **grz-pydantic-models:** read consent dateTimes as FHIR defines them ([#646](https://github.com/BfArM-MVH/grz-tools/issues/646)) ([c6c4f2f](https://github.com/BfArM-MVH/grz-tools/commit/c6c4f2ff051f581430653b2efcc4b156fe5ea33f)), closes [#645](https://github.com/BfArM-MVH/grz-tools/issues/645)
* **grzctl:** require grz-db 2.2 ([c6c4f2f](https://github.com/BfArM-MVH/grz-tools/commit/c6c4f2ff051f581430653b2efcc4b156fe5ea33f))
