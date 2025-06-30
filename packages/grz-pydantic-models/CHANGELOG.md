# Changelog

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
