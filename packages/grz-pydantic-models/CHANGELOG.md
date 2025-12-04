# Changelog

## [2.5.0](https://github.com/BfArM-MVH/grz-tools/compare/grz-pydantic-models-v2.4.0...grz-pydantic-models-v2.5.0) (2025-12-04)


### Features

* **grz-pydantic-models:** Add a thresholds model and implement check for oncomine panel submissions ([#469](https://github.com/BfArM-MVH/grz-tools/issues/469)) ([165a37c](https://github.com/BfArM-MVH/grz-tools/commit/165a37c697e1191c6c32d4162a32e531a18e2369)), closes [#462](https://github.com/BfArM-MVH/grz-tools/issues/462)
* **grz-pydantic-models:** add quarter date calculation functions ([d6a9618](https://github.com/BfArM-MVH/grz-tools/commit/d6a96182ecf7d827fe6d7dc12dd1efe3cf94c47a))
* **grz-pydantic-models:** lower long-read base quality thresholds to match BfArM ([#461](https://github.com/BfArM-MVH/grz-tools/issues/461)) ([d55800b](https://github.com/BfArM-MVH/grz-tools/commit/d55800b428e97b48d753b57fe59d1e66bb0aaeaf))
* **grzctl,grz-db:** add should-qc as db subcommand ([#473](https://github.com/BfArM-MVH/grz-tools/issues/473)) ([d6a9618](https://github.com/BfArM-MVH/grz-tools/commit/d6a96182ecf7d827fe6d7dc12dd1efe3cf94c47a))


### Bug Fixes

* **grz-pydantic-models:** Check for presence of necessary labData depending on the genomicStudySubtype ([#443](https://github.com/BfArM-MVH/grz-tools/issues/443)) ([6eac517](https://github.com/BfArM-MVH/grz-tools/commit/6eac517322efc3d0c40994f93d6d3c9fa7f55515))

## [2.4.0](https://github.com/BfArM-MVH/grz-tools/compare/grz-pydantic-models-v2.3.1...grz-pydantic-models-v2.4.0) (2025-10-23)


### Features

* **grz-common,grz-pydantic-models,grz-check:** use mean read length threshold ([#428](https://github.com/BfArM-MVH/grz-tools/issues/428)) ([18db996](https://github.com/BfArM-MVH/grz-tools/commit/18db99663f67b8883a038c61a765bcb1d2fb9edf))
* **grz-pydantic-models:** allow deprecated to be used with or without parentheses ([#418](https://github.com/BfArM-MVH/grz-tools/issues/418)) ([62b2ef3](https://github.com/BfArM-MVH/grz-tools/commit/62b2ef359ccc87cfea86a4f8ddb34e01e227e177))


### Bug Fixes

* **grz-pydantic-models:** ignore optional fields in required consent scope/category ([#429](https://github.com/BfArM-MVH/grz-tools/issues/429)) ([f9fad57](https://github.com/BfArM-MVH/grz-tools/commit/f9fad573c16056cd24c9c26db0fa6d10c2e4e755))

## [2.3.1](https://github.com/BfArM-MVH/grz-tools/compare/grz-pydantic-models-v2.3.0...grz-pydantic-models-v2.3.1) (2025-10-13)


### Bug Fixes

* **grz-pydantic-models:** add missing accepted versions ([#417](https://github.com/BfArM-MVH/grz-tools/issues/417)) ([ec85e8f](https://github.com/BfArM-MVH/grz-tools/commit/ec85e8f343af70cb6a74487bc8df1727459a3d48))
* **grz-pydantic-models:** don't generate submission ID from redacted TAN ([#412](https://github.com/BfArM-MVH/grz-tools/issues/412)) ([c925e53](https://github.com/BfArM-MVH/grz-tools/commit/c925e53d4c36e003bca62343913233c52af73b14))
* **grz-pydantic-models:** more conservative "deny" provisions for multiple codes ([#416](https://github.com/BfArM-MVH/grz-tools/issues/416)) ([dd36ab4](https://github.com/BfArM-MVH/grz-tools/commit/dd36ab42d4e6c232f510420eda00c94805af78f0))

## [2.3.0](https://github.com/BfArM-MVH/grz-tools/compare/grz-pydantic-models-v2.2.1...grz-pydantic-models-v2.3.0) (2025-10-07)


### Features

* **grz-pydantic-models:** accept submission metadata v1.3 in validation ([#395](https://github.com/BfArM-MVH/grz-tools/issues/395)) ([a3fc7ab](https://github.com/BfArM-MVH/grz-tools/commit/a3fc7ab27c8834a687fa1d7078aa8557e63e0d30))
* **grz-pydantic-models:** allow MV consent revocation on non-initial ([1179499](https://github.com/BfArM-MVH/grz-tools/commit/117949907151b612251ce5680d709d335f0e9427))
* **grz-pydantic-models:** require valid scope in metadata v1.3+ ([#396](https://github.com/BfArM-MVH/grz-tools/issues/396)) ([521b3a5](https://github.com/BfArM-MVH/grz-tools/commit/521b3a579464ff938c3abf2a0fb1e11e8f83d79a))
* **grz-pydantic-models:** support metadata schema v1.3 ([#378](https://github.com/BfArM-MVH/grz-tools/issues/378)) ([21a1ad5](https://github.com/BfArM-MVH/grz-tools/commit/21a1ad53eec40c554f6e1b3205620f6e0cb5033d))
* **grzctl:** add quarterly report export ([#376](https://github.com/BfArM-MVH/grz-tools/issues/376)) ([1179499](https://github.com/BfArM-MVH/grz-tools/commit/117949907151b612251ce5680d709d335f0e9427))


### Bug Fixes

* **grz-pydantic-models:** accept consent objects with no subprovisions ([#397](https://github.com/BfArM-MVH/grz-tools/issues/397)) ([916ff64](https://github.com/BfArM-MVH/grz-tools/commit/916ff644d54c51a24a0916afc23203b44bbe5150))

## [2.2.1](https://github.com/BfArM-MVH/grz-tools/compare/grz-pydantic-models-v2.2.0...grz-pydantic-models-v2.2.1) (2025-08-19)


### Bug Fixes

* **grz-pydantic-models:** require at least one DNA library for index donor ([#342](https://github.com/BfArM-MVH/grz-tools/issues/342)) ([14b3b93](https://github.com/BfArM-MVH/grz-tools/commit/14b3b934a4d221b004a8a781c51bc066f6798452))

## [2.2.0](https://github.com/BfArM-MVH/grz-tools/compare/grz-pydantic-models-v2.1.2...grz-pydantic-models-v2.2.0) (2025-07-31)


### Features

* **grzctl,grz-db,grz-common,grz-pydantic-models:** add columns, migration, and populate ([#306](https://github.com/BfArM-MVH/grz-tools/issues/306)) ([c158fa0](https://github.com/BfArM-MVH/grz-tools/commit/c158fa0cfe47ddacd66947dd57b814f43cfaefdc))

## [2.1.2](https://github.com/BfArM-MVH/grz-tools/compare/grz-pydantic-models-v2.1.1...grz-pydantic-models-v2.1.2) (2025-07-23)


### Bug Fixes

* **grz-pydantic-models:** warn if non-index donors do not pass thresholds ([#300](https://github.com/BfArM-MVH/grz-tools/issues/300)) ([00505a4](https://github.com/BfArM-MVH/grz-tools/commit/00505a486dad0062e7149135d63818a948c0e927))

## [2.1.1](https://github.com/BfArM-MVH/grz-tools/compare/grz-pydantic-models-v2.1.0...grz-pydantic-models-v2.1.1) (2025-07-18)


### Bug Fixes

* **grz-pydantic-models:** do not set a default schema URL ([#286](https://github.com/BfArM-MVH/grz-tools/issues/286)) ([4607dc2](https://github.com/BfArM-MVH/grz-tools/commit/4607dc2a0da0699594a3b5ff7ca219c1aca57638))
* **grz-pydantic-models:** downgrade VCF from required to recommended ([#289](https://github.com/BfArM-MVH/grz-tools/issues/289)) ([9e5d3d1](https://github.com/BfArM-MVH/grz-tools/commit/9e5d3d1a109eb9b422570c860c7d5272c372c177))

## [2.1.0](https://github.com/BfArM-MVH/grz-tools/compare/grz-pydantic-models-v2.0.3...grz-pydantic-models-v2.1.0) (2025-07-14)


### Features

* **grz-pydantic-models:** support test submission types and new schema versions ([#280](https://github.com/BfArM-MVH/grz-tools/issues/280)) ([e8f7701](https://github.com/BfArM-MVH/grz-tools/commit/e8f77013a31a4895d9a210eb348337e9725e8535))


### Bug Fixes

* **grz-pydantic-models:** Fix conditions for consent exemption ([#276](https://github.com/BfArM-MVH/grz-tools/issues/276)) ([bfe50a0](https://github.com/BfArM-MVH/grz-tools/commit/bfe50a040e0a4b7a2b0159f4c244cbe0b38ceeca))

## [2.0.3](https://github.com/BfArM-MVH/grz-tools/compare/grz-pydantic-models-v2.0.2...grz-pydantic-models-v2.0.3) (2025-07-09)


### Bug Fixes

* **grz-pydantic-models:** use correct key from thresholds definition ([#270](https://github.com/BfArM-MVH/grz-tools/issues/270)) ([485c504](https://github.com/BfArM-MVH/grz-tools/commit/485c504acd8648c3227182b5cbdf42195549554e))
* **grz-pydantic-models:** validate percentBasesAboveQualityThreshold against thresholds ([#268](https://github.com/BfArM-MVH/grz-tools/issues/268)) ([7df1f67](https://github.com/BfArM-MVH/grz-tools/commit/7df1f679760c399c8506453b7ef124b4dc142e60))

## [2.0.2](https://github.com/BfArM-MVH/grz-tools/compare/grz-pydantic-models-v2.0.1...grz-pydantic-models-v2.0.2) (2025-07-03)


### Bug Fixes

* **grz-pydantic-models,grz-cli,grzctl:** accept only metadata schema 1.1.9 ([#262](https://github.com/BfArM-MVH/grz-tools/issues/262)) ([f61bd2c](https://github.com/BfArM-MVH/grz-tools/commit/f61bd2c03e1a7ce0a667c7c9a7b467233d0835b3))

## [2.0.1](https://github.com/BfArM-MVH/grz-tools/compare/grz-pydantic-models-v2.0.0...grz-pydantic-models-v2.0.1) (2025-06-30)


### Bug Fixes

* **grz-pydantic-models:** add specific exemption to mvConsent ([#241](https://github.com/BfArM-MVH/grz-tools/issues/241)) ([779488d](https://github.com/BfArM-MVH/grz-tools/commit/779488d3c6f09b0c55c3faecfcd2698453874a64))

## [2.0.0](https://github.com/BfArM-MVH/grz-tools/compare/grz-pydantic-models-v1.5.0...grz-pydantic-models-v2.0.0) (2025-06-30)


### ⚠ BREAKING CHANGES

* **grz-pydantic-models:** drop deprecated functionality ([#220](https://github.com/BfArM-MVH/grz-tools/issues/220))

### Features

* **grz-pydantic-models,grzctl:** add optional research consent profile validation ([#165](https://github.com/BfArM-MVH/grz-tools/issues/165)) ([4a04dae](https://github.com/BfArM-MVH/grz-tools/commit/4a04daebf5936f0b398b2d7db03cf0f0f372970b))
* **grz-pydantic-models:** support submission metadata schema v1.1.9 ([#222](https://github.com/BfArM-MVH/grz-tools/issues/222)) ([5781a1b](https://github.com/BfArM-MVH/grz-tools/commit/5781a1b83a9e09a158a05862f107214c97d70994))


### Bug Fixes

* **grz-cli,grz-pydantic-models:** Disallow empty sequence data ([#218](https://github.com/BfArM-MVH/grz-tools/issues/218)) ([df28ab9](https://github.com/BfArM-MVH/grz-tools/commit/df28ab9dd78c97bdbbbcb68c4ff7a2208e049225))
* **grz-pydantic-models,grz-cli:** Check for duplicate checksums and file paths ([#182](https://github.com/BfArM-MVH/grz-tools/issues/182)) ([f01e705](https://github.com/BfArM-MVH/grz-tools/commit/f01e70595c232190a158906ba74ec180b4dcace9))
* **grz-pydantic-models,grz-common:** Allow symlinks ([#179](https://github.com/BfArM-MVH/grz-tools/issues/179)) ([43fcf7a](https://github.com/BfArM-MVH/grz-tools/commit/43fcf7ab1ae1a81aa79656073e764f310e5ed851))
* **grz-pydantic-models:** ensure filePath is normalized ([#217](https://github.com/BfArM-MVH/grz-tools/issues/217)) ([ffd8a9e](https://github.com/BfArM-MVH/grz-tools/commit/ffd8a9e1d6cbcd57ba5dc910a575ab5ba3ec651c))
* **grz-pydantic-models:** Ensure only one donor has relation 'index' ([#167](https://github.com/BfArM-MVH/grz-tools/issues/167)) ([9c48a1e](https://github.com/BfArM-MVH/grz-tools/commit/9c48a1ecdfcd10a8e15e9a55e79ea84be13c89c9))
* **grz-pydantic-models:** ensure unique donor pseudonyms within submission ([#181](https://github.com/BfArM-MVH/grz-tools/issues/181)) ([7f27037](https://github.com/BfArM-MVH/grz-tools/commit/7f27037c4fbc8ee8ccf1cb26ea15417a9dce70a4))
* **grz-pydantic-models:** ensure unique run IDs within a lab datum ([#231](https://github.com/BfArM-MVH/grz-tools/issues/231)) ([7f608fd](https://github.com/BfArM-MVH/grz-tools/commit/7f608fd7f43a8e596231a2bce1283cf29ef5a97c))
* **grz-pydantic-models:** prevent paired-end long read lab data ([#223](https://github.com/BfArM-MVH/grz-tools/issues/223)) ([e8979dc](https://github.com/BfArM-MVH/grz-tools/commit/e8979dc3fa83de229c1ccc091dcf35be957f781e))
* **grz-pydantic-models:** require valid file extensions for QC pipeline ([#158](https://github.com/BfArM-MVH/grz-tools/issues/158)) ([7fa69bd](https://github.com/BfArM-MVH/grz-tools/commit/7fa69bdcf6702a08c0b0409df37cec43d559f7ae))


### Miscellaneous Chores

* **grz-pydantic-models:** drop deprecated functionality ([#220](https://github.com/BfArM-MVH/grz-tools/issues/220)) ([a7a7e8e](https://github.com/BfArM-MVH/grz-tools/commit/a7a7e8e105c7eb2bb0d567b73bf4da76427fd4d3))

## [1.5.0](https://github.com/BfArM-MVH/grz-tools/compare/grz-pydantic-models-v1.4.0...grz-pydantic-models-v1.5.0) (2025-06-11)


### Features

* migrate to monorepo configuration ([36c7360](https://github.com/BfArM-MVH/grz-tools/commit/36c736044ce09473cc664b4471117465c5cab9a3))
