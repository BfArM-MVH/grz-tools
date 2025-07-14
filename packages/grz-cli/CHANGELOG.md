# Changelog

## [1.0.3](https://github.com/BfArM-MVH/grz-tools/compare/grz-cli-v1.0.2...grz-cli-v1.0.3) (2025-07-14)


### Bug Fixes

* **grzctl,grz-cli:** bump grz-pydantic-models for end-to-end tests ([#274](https://github.com/BfArM-MVH/grz-tools/issues/274)) ([a1ed791](https://github.com/BfArM-MVH/grz-tools/commit/a1ed791f1f9fce52d08f8e70fba12a674336d250))


### Documentation

* **grz-cli:** fix Conda update instructions ([#277](https://github.com/BfArM-MVH/grz-tools/issues/277)) ([fd6e054](https://github.com/BfArM-MVH/grz-tools/commit/fd6e054a8f8e80e1d3ab659ca1b9906199030c18))

## [1.0.2](https://github.com/BfArM-MVH/grz-tools/compare/grz-cli-v1.0.1...grz-cli-v1.0.2) (2025-07-04)


### Bug Fixes

* **grz-common,grzctl,grz-cli:** remove runtime dependency on type stubs ([#258](https://github.com/BfArM-MVH/grz-tools/issues/258)) ([a116499](https://github.com/BfArM-MVH/grz-tools/commit/a116499de19655ec9c4a43093c2c077dd10efbbc))
* **grz-db,grzctl:** address previously unchecked mypy type checks ([#247](https://github.com/BfArM-MVH/grz-tools/issues/247)) ([e51a65b](https://github.com/BfArM-MVH/grz-tools/commit/e51a65b090c891f44c6c4cc7199138d4cb15c07a))
* **grzctl,grz-cli:** bump deps to enforce 1.1.9 metadata submission ([#264](https://github.com/BfArM-MVH/grz-tools/issues/264)) ([1255260](https://github.com/BfArM-MVH/grz-tools/commit/1255260e4af25d342e1c17e803aa6f6152de69c7))
* **grzctl,grz-cli:** fix type check dependency by bumping grz-common ([#260](https://github.com/BfArM-MVH/grz-tools/issues/260)) ([92dce72](https://github.com/BfArM-MVH/grz-tools/commit/92dce723d8d2fbc7c11d03e2ebea98f7a0f4da19))

## [1.0.1](https://github.com/BfArM-MVH/grz-tools/compare/grz-cli-v1.0.0...grz-cli-v1.0.1) (2025-06-30)


### Bug Fixes

* **grz-cli,grzctl:** add specific exemption to mvConsent ([#243](https://github.com/BfArM-MVH/grz-tools/issues/243)) ([3ae4e55](https://github.com/BfArM-MVH/grz-tools/commit/3ae4e5513259933671146c40458e2c485a8fa612))

## [1.0.0](https://github.com/BfArM-MVH/grz-tools/compare/grz-cli-v0.7.0...grz-cli-v1.0.0) (2025-06-30)


### ⚠ BREAKING CHANGES

* **grz-cli,grzctl:** require GRZ and LE Id in config during validate ([#226](https://github.com/BfArM-MVH/grz-tools/issues/226))
* **grz-cli:** remove support for old configuration format ([#221](https://github.com/BfArM-MVH/grz-tools/issues/221))

### Features

* **grz-cli,grzctl:** require GRZ and LE Id in config during validate ([#226](https://github.com/BfArM-MVH/grz-tools/issues/226)) ([7043d9b](https://github.com/BfArM-MVH/grz-tools/commit/7043d9b3d66fcbd66bc102d9d0608467293ff7e1))


### Bug Fixes

* **grz-cli,grz-common:** Require click &gt;=8.2 ([#214](https://github.com/BfArM-MVH/grz-tools/issues/214)) ([bc6f839](https://github.com/BfArM-MVH/grz-tools/commit/bc6f839efa3a7b88025af66199b7eea06ac688ef))
* **grz-cli:** add missing validate parameter to submit command ([#235](https://github.com/BfArM-MVH/grz-tools/issues/235)) ([854e704](https://github.com/BfArM-MVH/grz-tools/commit/854e7046fc9bacf6a6aa71d8d6a0b4ea4b5d6d65))


### Miscellaneous Chores

* **grz-cli:** remove support for old configuration format ([#221](https://github.com/BfArM-MVH/grz-tools/issues/221)) ([49aeaae](https://github.com/BfArM-MVH/grz-tools/commit/49aeaae693a09a42505fe2606f1db8d85c1a4d21))

## [0.7.0](https://github.com/BfArM-MVH/grz-tools/compare/grz-cli-v0.6.1...grz-cli-v0.7.0) (2025-06-11)


### Features

* migrate to monorepo configuration ([36c7360](https://github.com/BfArM-MVH/grz-tools/commit/36c736044ce09473cc664b4471117465c5cab9a3))
