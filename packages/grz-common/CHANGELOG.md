# Changelog

## [1.0.3](https://github.com/BfArM-MVH/grz-tools/compare/grz-common-v1.0.2...grz-common-v1.0.3) (2025-07-14)


### Bug Fixes

* **grz-common:** reduce validation logging verbosity ([#273](https://github.com/BfArM-MVH/grz-tools/issues/273)) ([18d15d9](https://github.com/BfArM-MVH/grz-tools/commit/18d15d94543dcfb4ca3fa8918094ca1a52b8812c))

## [1.0.2](https://github.com/BfArM-MVH/grz-tools/compare/grz-common-v1.0.1...grz-common-v1.0.2) (2025-07-03)


### Bug Fixes

* **grz-common,grz-cli,grzctl:** fix encrypt caching ([#257](https://github.com/BfArM-MVH/grz-tools/issues/257)) ([3d86767](https://github.com/BfArM-MVH/grz-tools/commit/3d86767c77352e1a44807e312faac7604bd04de8))
* **grz-common,grzctl,grz-cli:** remove runtime dependency on type stubs ([#258](https://github.com/BfArM-MVH/grz-tools/issues/258)) ([a116499](https://github.com/BfArM-MVH/grz-tools/commit/a116499de19655ec9c4a43093c2c077dd10efbbc))

## [1.0.1](https://github.com/BfArM-MVH/grz-tools/compare/grz-common-v1.0.0...grz-common-v1.0.1) (2025-07-02)


### Bug Fixes

* **grz-cli,grzctl:** ignore missing log files if force enabled ([#249](https://github.com/BfArM-MVH/grz-tools/issues/249)) ([f29099f](https://github.com/BfArM-MVH/grz-tools/commit/f29099f147cbe5cdd1ad21eb5f3ef0e42d7385d2))
* **grz-db,grzctl:** address previously unchecked mypy type checks ([#247](https://github.com/BfArM-MVH/grz-tools/issues/247)) ([e51a65b](https://github.com/BfArM-MVH/grz-tools/commit/e51a65b090c891f44c6c4cc7199138d4cb15c07a))

## [1.0.0](https://github.com/BfArM-MVH/grz-tools/compare/grz-common-v0.1.0...grz-common-v1.0.0) (2025-06-30)


### ⚠ BREAKING CHANGES

* **grz-cli,grzctl:** require GRZ and LE Id in config during validate ([#226](https://github.com/BfArM-MVH/grz-tools/issues/226))

### Features

* **grz-cli,grzctl:** require GRZ and LE Id in config during validate ([#226](https://github.com/BfArM-MVH/grz-tools/issues/226)) ([7043d9b](https://github.com/BfArM-MVH/grz-tools/commit/7043d9b3d66fcbd66bc102d9d0608467293ff7e1))
* **grz-common:** check for existing files before encrypting ([#230](https://github.com/BfArM-MVH/grz-tools/issues/230)) ([28b84fd](https://github.com/BfArM-MVH/grz-tools/commit/28b84fd8a1133824c0ed624d494777d279f697eb))
* **grzctl,grz-common:** redact tanG and localCaseId before archiving ([#225](https://github.com/BfArM-MVH/grz-tools/issues/225)) ([b189b2c](https://github.com/BfArM-MVH/grz-tools/commit/b189b2ca94d59f2f640b07e0e6cc7e36df546049))
* **grzctl,grz-common:** use marker files while cleaning ([#228](https://github.com/BfArM-MVH/grz-tools/issues/228)) ([aacfaf9](https://github.com/BfArM-MVH/grz-tools/commit/aacfaf9a5da1c9d36835f679e522ef0376dde1d4))


### Bug Fixes

* **grz-cli,grz-common:** Require click &gt;=8.2 ([#214](https://github.com/BfArM-MVH/grz-tools/issues/214)) ([bc6f839](https://github.com/BfArM-MVH/grz-tools/commit/bc6f839efa3a7b88025af66199b7eea06ac688ef))
* **grz-cli,grzctl:** downgrade read length mismatch to warning ([#236](https://github.com/BfArM-MVH/grz-tools/issues/236)) ([bb3ebdb](https://github.com/BfArM-MVH/grz-tools/commit/bb3ebdb16b2baf4898e4683ed3c2c7eea9b07db2))
* **grz-common:** improve missing file error message ([#171](https://github.com/BfArM-MVH/grz-tools/issues/171)) ([4e110e7](https://github.com/BfArM-MVH/grz-tools/commit/4e110e7c96b387b2c4ae8390c400a5a6b004f2bb))
* **grz-common:** simplify read_multiple_json ([#180](https://github.com/BfArM-MVH/grz-tools/issues/180)) ([2abbb3c](https://github.com/BfArM-MVH/grz-tools/commit/2abbb3cb7d75d5d0a5b2fc85aaf10e83ad780793))
* **grz-pydantic-models,grz-common:** Allow symlinks ([#179](https://github.com/BfArM-MVH/grz-tools/issues/179)) ([43fcf7a](https://github.com/BfArM-MVH/grz-tools/commit/43fcf7ab1ae1a81aa79656073e764f310e5ed851))
* **grzctl:** rewrite download logic ([#183](https://github.com/BfArM-MVH/grz-tools/issues/183)) ([75894eb](https://github.com/BfArM-MVH/grz-tools/commit/75894ebbbbffd3125ae81a51927c1beff3b33990))

## 0.1.0 (2025-06-11)


### Features

* migrate to monorepo configuration ([36c7360](https://github.com/BfArM-MVH/grz-tools/commit/36c736044ce09473cc664b4471117465c5cab9a3))
