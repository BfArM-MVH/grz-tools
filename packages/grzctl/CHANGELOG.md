# Changelog

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
