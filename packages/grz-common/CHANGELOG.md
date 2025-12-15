# Changelog

## [1.6.1](https://github.com/BfArM-MVH/grz-tools/compare/grz-common-v1.6.0...grz-common-v1.6.1) (2025-12-14)


### Bug Fixes

* **grz-common,grzctl:** Revert submission date should be date of finished upload ([#460](https://github.com/BfArM-MVH/grz-tools/issues/460)) ([#478](https://github.com/BfArM-MVH/grz-tools/issues/478)) ([405793e](https://github.com/BfArM-MVH/grz-tools/commit/405793e1b40e5578d0f90d5d4d2daeb261d3da1f))

## [1.6.0](https://github.com/BfArM-MVH/grz-tools/compare/grz-common-v1.5.0...grz-common-v1.6.0) (2025-12-04)


### Features

* **grz-pydantic-models:** Add a thresholds model and implement check for oncomine panel submissions ([#469](https://github.com/BfArM-MVH/grz-tools/issues/469)) ([165a37c](https://github.com/BfArM-MVH/grz-tools/commit/165a37c697e1191c6c32d4162a32e531a18e2369)), closes [#462](https://github.com/BfArM-MVH/grz-tools/issues/462)


### Bug Fixes

* **grz-common,grzctl:** submission date should be date of finished upload ([#460](https://github.com/BfArM-MVH/grz-tools/issues/460)) ([907f194](https://github.com/BfArM-MVH/grz-tools/commit/907f194dda3b14ea2b69189febf3e490430f2a5d))

## [1.5.0](https://github.com/BfArM-MVH/grz-tools/compare/grz-common-v1.4.0...grz-common-v1.5.0) (2025-10-29)


### Features

* **grz-common:** upload version marker first ([#447](https://github.com/BfArM-MVH/grz-tools/issues/447)) ([156807f](https://github.com/BfArM-MVH/grz-tools/commit/156807f6cf00ef16d8f3223e83c4135b0591cdc7))
* **grzctl,grz-common:** show size in inbox listing ([#446](https://github.com/BfArM-MVH/grz-tools/issues/446)) ([3bee924](https://github.com/BfArM-MVH/grz-tools/commit/3bee92413f537fc4c087b0aca472548c557e3f78))

## [1.4.0](https://github.com/BfArM-MVH/grz-tools/compare/grz-common-v1.3.1...grz-common-v1.4.0) (2025-10-23)


### Features

* **grz-common,grz-pydantic-models,grz-check:** use mean read length threshold ([#428](https://github.com/BfArM-MVH/grz-tools/issues/428)) ([18db996](https://github.com/BfArM-MVH/grz-tools/commit/18db99663f67b8883a038c61a765bcb1d2fb9edf))


### Bug Fixes

* **grz-common:** meanReadLength param None vs 0 ([#435](https://github.com/BfArM-MVH/grz-tools/issues/435)) ([90538aa](https://github.com/BfArM-MVH/grz-tools/commit/90538aa6879b5f4e7c112cb7cb36ddb2ade3918c))
* **grz-common:** resolve KeyError with long-read submissions ([#432](https://github.com/BfArM-MVH/grz-tools/issues/432)) ([81bd0cb](https://github.com/BfArM-MVH/grz-tools/commit/81bd0cbf411282c69352eb86b3fc53258bfe9cc1))

## [1.3.1](https://github.com/BfArM-MVH/grz-tools/compare/grz-common-v1.3.0...grz-common-v1.3.1) (2025-10-13)


### Bug Fixes

* **grz-cli,grz-common:** Check encryption progress logs before upload ([#406](https://github.com/BfArM-MVH/grz-tools/issues/406)) ([401a20a](https://github.com/BfArM-MVH/grz-tools/commit/401a20aef1476eb940abc0c9aaf74e409215e55e))
* **grz-cli,grz-common:** check validation progress logs before encrypt ([#411](https://github.com/BfArM-MVH/grz-tools/issues/411)) ([a33b342](https://github.com/BfArM-MVH/grz-tools/commit/a33b342da6f32a59c57acfa067584ac7798f9764))
* **grz-pydantic-models:** don't generate submission ID from redacted TAN ([#412](https://github.com/BfArM-MVH/grz-tools/issues/412)) ([c925e53](https://github.com/BfArM-MVH/grz-tools/commit/c925e53d4c36e003bca62343913233c52af73b14))

## [1.3.0](https://github.com/BfArM-MVH/grz-tools/compare/grz-common-v1.2.1...grz-common-v1.3.0) (2025-08-19)


### Features

* **grz-common:** add extra warnings for BAM files ([#356](https://github.com/BfArM-MVH/grz-tools/issues/356)) ([6c3b62a](https://github.com/BfArM-MVH/grz-tools/commit/6c3b62a6ada98e7a7ad4bb15b6ad86cf7c27c5c1))
* **grz-common:** Add proxy options to config ([#335](https://github.com/BfArM-MVH/grz-tools/issues/335)) ([a9acc4e](https://github.com/BfArM-MVH/grz-tools/commit/a9acc4ee487de19cc29965aff1d1b10a32f174f7))

## [1.2.1](https://github.com/BfArM-MVH/grz-tools/compare/grz-common-v1.2.0...grz-common-v1.2.1) (2025-08-05)


### Bug Fixes

* **grz-cli,grz-common,grzctl:** fix logging and migrate to grz-common ([#319](https://github.com/BfArM-MVH/grz-tools/issues/319)) ([51ada07](https://github.com/BfArM-MVH/grz-tools/commit/51ada073a2af93ba1a1c48f069b4546ce9bd2975))

## [1.2.0](https://github.com/BfArM-MVH/grz-tools/compare/grz-common-v1.1.1...grz-common-v1.2.0) (2025-07-31)


### Features

* **grzctl,grz-db,grz-common,grz-pydantic-models:** add columns, migration, and populate ([#306](https://github.com/BfArM-MVH/grz-tools/issues/306)) ([c158fa0](https://github.com/BfArM-MVH/grz-tools/commit/c158fa0cfe47ddacd66947dd57b814f43cfaefdc))


### Bug Fixes

* **grz-common:** bump grz-pydanic-models version ([#316](https://github.com/BfArM-MVH/grz-tools/issues/316)) ([b4cd8e2](https://github.com/BfArM-MVH/grz-tools/commit/b4cd8e2925a24e7822ede3ddbfa1def4dccf8b87))

## [1.1.1](https://github.com/BfArM-MVH/grz-tools/compare/grz-common-v1.1.0...grz-common-v1.1.1) (2025-07-23)


### Bug Fixes

* **grz-common:** correctly fallback without grz-check ([#298](https://github.com/BfArM-MVH/grz-tools/issues/298)) ([1ca0b20](https://github.com/BfArM-MVH/grz-tools/commit/1ca0b20a7f2b57f5e49144956fb86d9de1b4301c))

## [1.1.0](https://github.com/BfArM-MVH/grz-tools/compare/grz-common-v1.0.3...grz-common-v1.1.0) (2025-07-21)


### Features

* **grz-cli,grzctl,grz-check:** Add grz-check ([#240](https://github.com/BfArM-MVH/grz-tools/issues/240)) ([57048f6](https://github.com/BfArM-MVH/grz-tools/commit/57048f66888cb566887e627a2b973c3f8b1b83c5))

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
