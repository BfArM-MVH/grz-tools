# Changelog

## [0.5.0](https://github.com/BfArM-MVH/grz-tools/compare/grzctl-v0.4.0...grzctl-v0.5.0) (2025-08-19)


### Features

* **grz-cli,grzctl:** bump grz dependencies to latest ([#359](https://github.com/BfArM-MVH/grz-tools/issues/359)) ([b31c88b](https://github.com/BfArM-MVH/grz-tools/commit/b31c88bfa32cf257e6db9b2bb80d302493dbd469))
* **grz-cli:** also report grz library versions ([#355](https://github.com/BfArM-MVH/grz-tools/issues/355)) ([58d047e](https://github.com/BfArM-MVH/grz-tools/commit/58d047ec3492df172067312203da49c8af83a1a8))
* **grz-db:** add limit parameter to database list_submissions ([e2eebda](https://github.com/BfArM-MVH/grz-tools/commit/e2eebdaaaa524cfeacb97f9717ba85bd74b2c8a6))
* **grzctl:** add configurable display limit for db list command ([#344](https://github.com/BfArM-MVH/grz-tools/issues/344)) ([e2eebda](https://github.com/BfArM-MVH/grz-tools/commit/e2eebdaaaa524cfeacb97f9717ba85bd74b2c8a6))
* **grzctl:** add configurable display limit to list command ([9279364](https://github.com/BfArM-MVH/grz-tools/commit/927936449f3d9e1707e631f343222b53b33b9862))
* **grzctl:** also report grz library versions ([58d047e](https://github.com/BfArM-MVH/grz-tools/commit/58d047ec3492df172067312203da49c8af83a1a8))
* **grzctl:** confirm before updating submission from error state ([#357](https://github.com/BfArM-MVH/grz-tools/issues/357)) ([25e6cb6](https://github.com/BfArM-MVH/grz-tools/commit/25e6cb62130cf926a9c77d5232bc39d3ecb91c66))
* **grzctl:** optionally add database state to list output ([#341](https://github.com/BfArM-MVH/grz-tools/issues/341)) ([9279364](https://github.com/BfArM-MVH/grz-tools/commit/927936449f3d9e1707e631f343222b53b33b9862))
* **grzctl:** replace oldest upload with duration in inbox list ([#345](https://github.com/BfArM-MVH/grz-tools/issues/345)) ([e16f796](https://github.com/BfArM-MVH/grz-tools/commit/e16f796251db490a9e08b79ac01b1e110e3a318b))


### Bug Fixes

* **grz-db:** allow empty author private key passphrases ([25e6cb6](https://github.com/BfArM-MVH/grz-tools/commit/25e6cb62130cf926a9c77d5232bc39d3ecb91c66))
* **grzctl:** `pruefbericht --dry-run` does not require ([4183595](https://github.com/BfArM-MVH/grz-tools/commit/4183595dabf9324fe81a2767b93c5be03b674d60))
* **grzctl:** `pruefbericht --dry-run` does not require config/credentials ([#358](https://github.com/BfArM-MVH/grz-tools/issues/358)) ([4183595](https://github.com/BfArM-MVH/grz-tools/commit/4183595dabf9324fe81a2767b93c5be03b674d60)), closes [#176](https://github.com/BfArM-MVH/grz-tools/issues/176)
* **grzctl:** display database state column if table empty but config ([e16f796](https://github.com/BfArM-MVH/grz-tools/commit/e16f796251db490a9e08b79ac01b1e110e3a318b))

## [0.4.0](https://github.com/BfArM-MVH/grz-tools/compare/grzctl-v0.3.0...grzctl-v0.4.0) (2025-08-05)


### Features

* **grzctl:** add reporting for processed submissions ([#320](https://github.com/BfArM-MVH/grz-tools/issues/320)) ([d44aead](https://github.com/BfArM-MVH/grz-tools/commit/d44aeade809e39693360b577e5482873ae975709))


### Bug Fixes

* **grz-cli,grz-common,grzctl:** fix logging and migrate to grz-common ([#319](https://github.com/BfArM-MVH/grz-tools/issues/319)) ([51ada07](https://github.com/BfArM-MVH/grz-tools/commit/51ada073a2af93ba1a1c48f069b4546ce9bd2975))

## [0.3.0](https://github.com/BfArM-MVH/grz-tools/compare/grzctl-v0.2.6...grzctl-v0.3.0) (2025-07-31)


### Features

* **grzctl,grz-db,grz-common,grz-pydantic-models:** add columns, migration, and populate ([#306](https://github.com/BfArM-MVH/grz-tools/issues/306)) ([c158fa0](https://github.com/BfArM-MVH/grz-tools/commit/c158fa0cfe47ddacd66947dd57b814f43cfaefdc))


### Bug Fixes

* **grzctl:** bump dependencies for migrations ([#317](https://github.com/BfArM-MVH/grz-tools/issues/317)) ([3f2e529](https://github.com/BfArM-MVH/grz-tools/commit/3f2e52976bcedfb7e355b58972e7860b586243af))

## [0.2.6](https://github.com/BfArM-MVH/grz-tools/compare/grzctl-v0.2.5...grzctl-v0.2.6) (2025-07-23)


### Bug Fixes

* **grzctl,grz-cli:** bump dependencies ([#303](https://github.com/BfArM-MVH/grz-tools/issues/303)) ([cf12d35](https://github.com/BfArM-MVH/grz-tools/commit/cf12d35a7a20dcb5494a3576ccc06c393f763367))

## [0.2.5](https://github.com/BfArM-MVH/grz-tools/compare/grzctl-v0.2.4...grzctl-v0.2.5) (2025-07-21)


### Bug Fixes

* **grz-cli,grzctl:** bump deps to make VCFs optional ([#291](https://github.com/BfArM-MVH/grz-tools/issues/291)) ([fe4ac2f](https://github.com/BfArM-MVH/grz-tools/commit/fe4ac2f8230b804ad9ec2c6d2102207ab97b0365))

## [0.2.4](https://github.com/BfArM-MVH/grz-tools/compare/grzctl-v0.2.3...grzctl-v0.2.4) (2025-07-14)


### Bug Fixes

* **grzctl,grz-cli:** bump grz-pydantic-models for end-to-end tests ([#274](https://github.com/BfArM-MVH/grz-tools/issues/274)) ([a1ed791](https://github.com/BfArM-MVH/grz-tools/commit/a1ed791f1f9fce52d08f8e70fba12a674336d250))

## [0.2.3](https://github.com/BfArM-MVH/grz-tools/compare/grzctl-v0.2.2...grzctl-v0.2.3) (2025-07-04)


### Bug Fixes

* **grz-common,grzctl,grz-cli:** remove runtime dependency on type stubs ([#258](https://github.com/BfArM-MVH/grz-tools/issues/258)) ([a116499](https://github.com/BfArM-MVH/grz-tools/commit/a116499de19655ec9c4a43093c2c077dd10efbbc))
* **grzctl,grz-cli:** bump deps to enforce 1.1.9 metadata submission ([#264](https://github.com/BfArM-MVH/grz-tools/issues/264)) ([1255260](https://github.com/BfArM-MVH/grz-tools/commit/1255260e4af25d342e1c17e803aa6f6152de69c7))
* **grzctl,grz-cli:** fix type check dependency by bumping grz-common ([#260](https://github.com/BfArM-MVH/grz-tools/issues/260)) ([92dce72](https://github.com/BfArM-MVH/grz-tools/commit/92dce723d8d2fbc7c11d03e2ebea98f7a0f4da19))

## [0.2.2](https://github.com/BfArM-MVH/grz-tools/compare/grzctl-v0.2.1...grzctl-v0.2.2) (2025-07-02)


### Bug Fixes

* **grz-db,grzctl:** address previously unchecked mypy type checks ([#247](https://github.com/BfArM-MVH/grz-tools/issues/247)) ([e51a65b](https://github.com/BfArM-MVH/grz-tools/commit/e51a65b090c891f44c6c4cc7199138d4cb15c07a))
* **grz-db:** undefined verifying_key_comment if no latest_state_obj ([#246](https://github.com/BfArM-MVH/grz-tools/issues/246)) ([dc793aa](https://github.com/BfArM-MVH/grz-tools/commit/dc793aaa4be33ff2a55dd2017869dc3bfea9f22d))

## [0.2.1](https://github.com/BfArM-MVH/grz-tools/compare/grzctl-v0.2.0...grzctl-v0.2.1) (2025-06-30)


### Bug Fixes

* **grz-cli,grzctl:** add specific exemption to mvConsent ([#243](https://github.com/BfArM-MVH/grz-tools/issues/243)) ([3ae4e55](https://github.com/BfArM-MVH/grz-tools/commit/3ae4e5513259933671146c40458e2c485a8fa612))

## [0.2.0](https://github.com/BfArM-MVH/grz-tools/compare/grzctl-v0.1.0...grzctl-v0.2.0) (2025-06-30)


### Features

* **grz-pydantic-models,grzctl:** add optional research consent profile validation ([#165](https://github.com/BfArM-MVH/grz-tools/issues/165)) ([4a04dae](https://github.com/BfArM-MVH/grz-tools/commit/4a04daebf5936f0b398b2d7db03cf0f0f372970b))
* **grzctl,grz-common:** use marker files while cleaning ([#228](https://github.com/BfArM-MVH/grz-tools/issues/228)) ([aacfaf9](https://github.com/BfArM-MVH/grz-tools/commit/aacfaf9a5da1c9d36835f679e522ef0376dde1d4))
* **grzctl,grz-db:** Add support for change requests ([#151](https://github.com/BfArM-MVH/grz-tools/issues/151)) ([2f28d69](https://github.com/BfArM-MVH/grz-tools/commit/2f28d691b72da2d904391680ff72b1f9a3a22254))
* **grzctl:** Add --dry-run option to print pruefbericht, then exit ([#169](https://github.com/BfArM-MVH/grz-tools/issues/169)) ([1cea450](https://github.com/BfArM-MVH/grz-tools/commit/1cea4500dd563fe46531e535471afcb8b3b1bb8e))


### Bug Fixes

* **grz-cli,grz-common:** Require click &gt;=8.2 ([#214](https://github.com/BfArM-MVH/grz-tools/issues/214)) ([bc6f839](https://github.com/BfArM-MVH/grz-tools/commit/bc6f839efa3a7b88025af66199b7eea06ac688ef))
* **grz-db,grzctl:** constrain author names slightly ([#205](https://github.com/BfArM-MVH/grz-tools/issues/205)) ([2d46420](https://github.com/BfArM-MVH/grz-tools/commit/2d464204fb8d07773d04d31b6fa93208e4181f22))
* **grzctl,grz-db:** Add `db submission modify` to allow setting tanG/pseudonym ([#198](https://github.com/BfArM-MVH/grz-tools/issues/198)) ([b6275c3](https://github.com/BfArM-MVH/grz-tools/commit/b6275c38b134e6d334dc158c9c98631e62750b68))
* **grzctl:** add missing grz-cli dependency ([#154](https://github.com/BfArM-MVH/grz-tools/issues/154)) ([57b77ad](https://github.com/BfArM-MVH/grz-tools/commit/57b77adb358ac1c5befba6df3dbc4297908fb953))
* **grzctl:** disallow empty submission id for `clean` ([#186](https://github.com/BfArM-MVH/grz-tools/issues/186)) ([c9556a4](https://github.com/BfArM-MVH/grz-tools/commit/c9556a4d3d80dbf36d2e5d7feed573df248711b4))
* **grzctl:** submit proper library type in Prüfbericht ([#219](https://github.com/BfArM-MVH/grz-tools/issues/219)) ([45ebc14](https://github.com/BfArM-MVH/grz-tools/commit/45ebc14e408558f06fd92055f74efc092002174a))


### Documentation

* **grzctl:** add instructions for running dev version ([#213](https://github.com/BfArM-MVH/grz-tools/issues/213)) ([2557b24](https://github.com/BfArM-MVH/grz-tools/commit/2557b24885a04c85ef1156b56471d98c933ff81d))

## 0.1.0 (2025-06-11)


### Features

* migrate to monorepo configuration ([36c7360](https://github.com/BfArM-MVH/grz-tools/commit/36c736044ce09473cc664b4471117465c5cab9a3))
