# Changelog

## [2.1.3](https://github.com/BfArM-MVH/grz-tools/compare/grz-db-v2.1.2...grz-db-v2.1.3) (2026-07-21)


### Bug Fixes

* **grz-common,grz-db:** map duplicate tanG and upload errors correctly   ([#617](https://github.com/BfArM-MVH/grz-tools/issues/617)) ([7c6fe27](https://github.com/BfArM-MVH/grz-tools/commit/7c6fe271a15011a6925a2f2a8a1167f832d78f9c))
* **grz-db:** datetime handling in signature evaluation ([#626](https://github.com/BfArM-MVH/grz-tools/issues/626)) ([b7894c9](https://github.com/BfArM-MVH/grz-tools/commit/b7894c9bee81acff0389bdef9e7414bb1ba795b9))

## [2.1.2](https://github.com/BfArM-MVH/grz-tools/compare/grz-db-v2.1.1...grz-db-v2.1.2) (2026-06-08)


### Bug Fixes

* **grz-db:** correct tan_g population by raising sqlmodel floor ([#594](https://github.com/BfArM-MVH/grz-tools/issues/594)) ([abaa69f](https://github.com/BfArM-MVH/grz-tools/commit/abaa69f85ce84c0e3f69feaac0223b2622d16916))

## [2.1.1](https://github.com/BfArM-MVH/grz-tools/compare/grz-db-v2.1.0...grz-db-v2.1.1) (2026-05-26)


### Bug Fixes

* **grz-pydantic-models:** force timezone during consent_by_code ([98571dd](https://github.com/BfArM-MVH/grz-tools/commit/98571dd409e98dedb29d74360d21ee726f3f5617))
* **grzctl:** fix backfill datetime vs date issues ([98571dd](https://github.com/BfArM-MVH/grz-tools/commit/98571dd409e98dedb29d74360d21ee726f3f5617))
* **grzctl:** fix backfill datetime vs date issues ([#583](https://github.com/BfArM-MVH/grz-tools/issues/583)) ([98571dd](https://github.com/BfArM-MVH/grz-tools/commit/98571dd409e98dedb29d74360d21ee726f3f5617))

## [2.1.0](https://github.com/BfArM-MVH/grz-tools/compare/grz-db-v2.0.0...grz-db-v2.1.0) (2026-05-19)


### Features

* **grz-db,grz-common:** clean up SubmissionDb.populate API surface ([a0fc796](https://github.com/BfArM-MVH/grz-tools/commit/a0fc7961e30c7689d469ecf4fcb3eec65f4a8398))


### Bug Fixes

* **grz-common,grz-db,grzctl:** move populate logic into SubmissionDb ([a0fc796](https://github.com/BfArM-MVH/grz-tools/commit/a0fc7961e30c7689d469ecf4fcb3eec65f4a8398))
* **grz-pydantic-models,grz-db:** align inconsistent version pins ([a0fc796](https://github.com/BfArM-MVH/grz-tools/commit/a0fc7961e30c7689d469ecf4fcb3eec65f4a8398))
* **grzctl:** correct stale grz-cli hint in submission-not-found errors ([a0fc796](https://github.com/BfArM-MVH/grz-tools/commit/a0fc7961e30c7689d469ecf4fcb3eec65f4a8398))

## [2.0.0](https://github.com/BfArM-MVH/grz-tools/compare/grz-db-v1.2.0...grz-db-v2.0.0) (2026-05-11)


### ⚠ BREAKING CHANGES

* **grzctl:** Automatically update submission state upon subcommand invocation ([#501](https://github.com/BfArM-MVH/grz-tools/issues/501))

### Features

* **grz-cli, grz-db:** require & record QC workflow version; log grzctl runtime version (closes [#532](https://github.com/BfArM-MVH/grz-tools/issues/532)) ([#561](https://github.com/BfArM-MVH/grz-tools/issues/561)) ([b0929cf](https://github.com/BfArM-MVH/grz-tools/commit/b0929cfd8b77fbc8d6f9e940ff8d2834b1c3271e))
* **grz-db,grzctl:** add selected_for_qc submission column to submission database ([#531](https://github.com/BfArM-MVH/grz-tools/issues/531)) ([a753bd2](https://github.com/BfArM-MVH/grz-tools/commit/a753bd2f44a73fc113e10a736dd8c09f7c035b07))
* **grz-tools:** add structured failure reason tracking for submissions ([#544](https://github.com/BfArM-MVH/grz-tools/issues/544)) ([65a2143](https://github.com/BfArM-MVH/grz-tools/commit/65a214389e91be00e99fb812293bfe40cca45cf0))
* **grzctl, grz-db, grz-pydantic-models:** Improved populate ([#547](https://github.com/BfArM-MVH/grz-tools/issues/547)) ([8cbf232](https://github.com/BfArM-MVH/grz-tools/commit/8cbf2323747a3a2c35623b80f7ad1ef420217f07))
* **grzctl:** add state-based filtering to grzctl db list (latest by default) ([#523](https://github.com/BfArM-MVH/grz-tools/issues/523)) ([9a715f4](https://github.com/BfArM-MVH/grz-tools/commit/9a715f40ec80d149fe5ab2fbda383da9671ff332)), closes [#504](https://github.com/BfArM-MVH/grz-tools/issues/504)
* **grzctl:** Automatically update submission state upon subcommand invocation ([#501](https://github.com/BfArM-MVH/grz-tools/issues/501)) ([8953102](https://github.com/BfArM-MVH/grz-tools/commit/895310211d36796d995c37f1929a12458195d4c7))
* **grzctl:** Introduce "backfill" command to re-read metadata from the archive ([#559](https://github.com/BfArM-MVH/grz-tools/issues/559)) ([9e83ad3](https://github.com/BfArM-MVH/grz-tools/commit/9e83ad31dbb9dd44f18d1f0ed1dbef82f4f1f63e))


### Bug Fixes

* **repo:** Update dependencies ([#498](https://github.com/BfArM-MVH/grz-tools/issues/498)) ([368dfdb](https://github.com/BfArM-MVH/grz-tools/commit/368dfdbaa703f17f0c290ea051be30f9be4bebf3))

## [1.2.0](https://github.com/BfArM-MVH/grz-tools/compare/grz-db-v1.1.0...grz-db-v1.2.0) (2025-12-04)


### Features

* **grz-db,grzctl:** add PostgreSQL support ([#459](https://github.com/BfArM-MVH/grz-tools/issues/459)) ([9c0f941](https://github.com/BfArM-MVH/grz-tools/commit/9c0f941d156d19305d1603d65793f1a7dfda4756))
* **grz-pydantic-models:** add quarter date calculation functions ([d6a9618](https://github.com/BfArM-MVH/grz-tools/commit/d6a96182ecf7d827fe6d7dc12dd1efe3cf94c47a))
* **grzctl,grz-db:** add should-qc as db subcommand ([#473](https://github.com/BfArM-MVH/grz-tools/issues/473)) ([d6a9618](https://github.com/BfArM-MVH/grz-tools/commit/d6a96182ecf7d827fe6d7dc12dd1efe3cf94c47a))


### Bug Fixes

* **grz-db:** fix relation SQL Enum values ([#465](https://github.com/BfArM-MVH/grz-tools/issues/465)) ([8bb7330](https://github.com/BfArM-MVH/grz-tools/commit/8bb733001b50df6d645d21d252dad91917afebf1))

## [1.1.0](https://github.com/BfArM-MVH/grz-tools/compare/grz-db-v1.0.1...grz-db-v1.1.0) (2025-10-23)


### Features

* **grzctl,grz-cli:** support meanReadLength ([#437](https://github.com/BfArM-MVH/grz-tools/issues/437)) ([b86b843](https://github.com/BfArM-MVH/grz-tools/commit/b86b84313758d6fa16b1ee74af4834ba3e2ec914))


### Bug Fixes

* **grz-db:** bump grz-pydantic-models version ([b86b843](https://github.com/BfArM-MVH/grz-tools/commit/b86b84313758d6fa16b1ee74af4834ba3e2ec914))

## [1.0.1](https://github.com/BfArM-MVH/grz-tools/compare/grz-db-v1.0.0...grz-db-v1.0.1) (2025-10-13)


### Bug Fixes

* **grz-db,grzctl:** properly repopulate donors ([#413](https://github.com/BfArM-MVH/grz-tools/issues/413)) ([cb7b1bd](https://github.com/BfArM-MVH/grz-tools/commit/cb7b1bdebcfaec2e5581eb8a2e93ab57397b242e))

## [1.0.0](https://github.com/BfArM-MVH/grz-tools/compare/grz-db-v0.6.0...grz-db-v1.0.0) (2025-10-07)


### Features

* **grz-db:** bump required dependencies ([#401](https://github.com/BfArM-MVH/grz-tools/issues/401)) ([f62a6e1](https://github.com/BfArM-MVH/grz-tools/commit/f62a6e1982f7cd43210d9abf1856f7a46607092c))
* **grz-pydantic-models:** allow MV consent revocation on non-initial ([1179499](https://github.com/BfArM-MVH/grz-tools/commit/117949907151b612251ce5680d709d335f0e9427))
* **grzctl:** add quarterly report export ([#376](https://github.com/BfArM-MVH/grz-tools/issues/376)) ([1179499](https://github.com/BfArM-MVH/grz-tools/commit/117949907151b612251ce5680d709d335f0e9427))

## [0.6.0](https://github.com/BfArM-MVH/grz-tools/compare/grz-db-v0.5.0...grz-db-v0.6.0) (2025-08-27)


### Features

* **grz-db,grzctl:** sort submission db list by latest state with fallbacks ([#370](https://github.com/BfArM-MVH/grz-tools/issues/370)) ([fdc521b](https://github.com/BfArM-MVH/grz-tools/commit/fdc521bcc28af3c036aea7fa89837fa078eec25f))

## [0.5.0](https://github.com/BfArM-MVH/grz-tools/compare/grz-db-v0.4.0...grz-db-v0.5.0) (2025-08-19)


### Features

* **grz-db:** add limit parameter to database list_submissions ([e2eebda](https://github.com/BfArM-MVH/grz-tools/commit/e2eebdaaaa524cfeacb97f9717ba85bd74b2c8a6))
* **grzctl:** add configurable display limit for db list command ([#344](https://github.com/BfArM-MVH/grz-tools/issues/344)) ([e2eebda](https://github.com/BfArM-MVH/grz-tools/commit/e2eebdaaaa524cfeacb97f9717ba85bd74b2c8a6))
* **grzctl:** confirm before updating submission from error state ([#357](https://github.com/BfArM-MVH/grz-tools/issues/357)) ([25e6cb6](https://github.com/BfArM-MVH/grz-tools/commit/25e6cb62130cf926a9c77d5232bc39d3ecb91c66))


### Bug Fixes

* **grz-db:** allow empty author private key passphrases ([25e6cb6](https://github.com/BfArM-MVH/grz-tools/commit/25e6cb62130cf926a9c77d5232bc39d3ecb91c66))

## [0.4.0](https://github.com/BfArM-MVH/grz-tools/compare/grz-db-v0.3.0...grz-db-v0.4.0) (2025-08-05)


### Features

* **grzctl:** add reporting for processed submissions ([#320](https://github.com/BfArM-MVH/grz-tools/issues/320)) ([d44aead](https://github.com/BfArM-MVH/grz-tools/commit/d44aeade809e39693360b577e5482873ae975709))


### Bug Fixes

* **grz-db:** add StringConstraints to SubmissionBase.id ([#323](https://github.com/BfArM-MVH/grz-tools/issues/323)) ([0ac80fb](https://github.com/BfArM-MVH/grz-tools/commit/0ac80fbb4e68957bb9b59a395c90bc2bdf67e02d))

## [0.3.0](https://github.com/BfArM-MVH/grz-tools/compare/grz-db-v0.2.1...grz-db-v0.3.0) (2025-07-31)


### Features

* **grzctl,grz-db,grz-common,grz-pydantic-models:** add columns, migration, and populate ([#306](https://github.com/BfArM-MVH/grz-tools/issues/306)) ([c158fa0](https://github.com/BfArM-MVH/grz-tools/commit/c158fa0cfe47ddacd66947dd57b814f43cfaefdc))

## [0.2.1](https://github.com/BfArM-MVH/grz-tools/compare/grz-db-v0.2.0...grz-db-v0.2.1) (2025-07-02)


### Bug Fixes

* **grz-db,grzctl:** address previously unchecked mypy type checks ([#247](https://github.com/BfArM-MVH/grz-tools/issues/247)) ([e51a65b](https://github.com/BfArM-MVH/grz-tools/commit/e51a65b090c891f44c6c4cc7199138d4cb15c07a))

## [0.2.0](https://github.com/BfArM-MVH/grz-tools/compare/grz-db-v0.1.0...grz-db-v0.2.0) (2025-06-30)


### Features

* **grzctl,grz-db:** Add support for change requests ([#151](https://github.com/BfArM-MVH/grz-tools/issues/151)) ([2f28d69](https://github.com/BfArM-MVH/grz-tools/commit/2f28d691b72da2d904391680ff72b1f9a3a22254))


### Bug Fixes

* **grz-db:** don't print the duplicate tanG ([#229](https://github.com/BfArM-MVH/grz-tools/issues/229)) ([718a2be](https://github.com/BfArM-MVH/grz-tools/commit/718a2be52d959be44449f6b46143be62728c2631))
* **grz-db:** Ensure author's key is ed25519 ([#204](https://github.com/BfArM-MVH/grz-tools/issues/204)) ([0f7eba2](https://github.com/BfArM-MVH/grz-tools/commit/0f7eba2652c67f3c4ddb507f7d4e197dc0c086ec))
* **grzctl,grz-db:** Add `db submission modify` to allow setting tanG/pseudonym ([#198](https://github.com/BfArM-MVH/grz-tools/issues/198)) ([b6275c3](https://github.com/BfArM-MVH/grz-tools/commit/b6275c38b134e6d334dc158c9c98631e62750b68))
* **grzctl:** improve error message on incorrect passphrase for private key ([#206](https://github.com/BfArM-MVH/grz-tools/issues/206)) ([8b73036](https://github.com/BfArM-MVH/grz-tools/commit/8b7303643b96b87bf9b095e135633fc3db3a7c7e))

## 0.1.0 (2025-06-11)


### Features

* migrate to monorepo configuration ([36c7360](https://github.com/BfArM-MVH/grz-tools/commit/36c736044ce09473cc664b4471117465c5cab9a3))
