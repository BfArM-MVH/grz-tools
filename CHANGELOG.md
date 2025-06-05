# Changelog

## [0.6.1](https://github.com/BfArM-MVH/grz-cli/compare/v0.6.0...v0.6.1) (2025-06-05)


### Bug Fixes

* threads parameter was ignored in upload and download ([#118](https://github.com/BfArM-MVH/grz-cli/issues/118)) ([22692d6](https://github.com/BfArM-MVH/grz-cli/commit/22692d6be756636b5ab14d81a4a51cfb628d0cd1))

## [0.6.0](https://github.com/BfArM-MVH/grz-cli/compare/v0.5.0...v0.6.0) (2025-05-27)


### Features

* add `consent` command ([#97](https://github.com/BfArM-MVH/grz-cli/issues/97)) ([1074adf](https://github.com/BfArM-MVH/grz-cli/commit/1074adf65f00025b288a606d1a905137b93ce52f))
* add archive command ([#104](https://github.com/BfArM-MVH/grz-cli/issues/104)) ([a68b976](https://github.com/BfArM-MVH/grz-cli/commit/a68b976631b35be2ced260728af57846d62aecc1))
* add command to submit Prüfbericht ([#98](https://github.com/BfArM-MVH/grz-cli/issues/98)) ([e7841aa](https://github.com/BfArM-MVH/grz-cli/commit/e7841aab3034eaf27144a9db3e27821bc3f2106d))
* add force flag to ignore cached results ([#95](https://github.com/BfArM-MVH/grz-cli/issues/95)) ([fbdff96](https://github.com/BfArM-MVH/grz-cli/commit/fbdff961685a45de80aea5b9f430f57eded42ffc))
* cap filename length to prevent upload errors ([#96](https://github.com/BfArM-MVH/grz-cli/issues/96)) ([dd0be79](https://github.com/BfArM-MVH/grz-cli/commit/dd0be79f4836ba7b0d7a75679bbd643dd9c4cef8))
* declare typing support for entire package ([f05a4c7](https://github.com/BfArM-MVH/grz-cli/commit/f05a4c704fbf8f3395fa19d1f3b66c8f8a2669e9))
* implement `grz-cli clean` cmd ([#90](https://github.com/BfArM-MVH/grz-cli/issues/90)) ([4397fd0](https://github.com/BfArM-MVH/grz-cli/commit/4397fd09be613a1f7693e128411c4321638dca71))


### Bug Fixes

* allow listing empty buckets ([#103](https://github.com/BfArM-MVH/grz-cli/issues/103)) ([4f2e9ee](https://github.com/BfArM-MVH/grz-cli/commit/4f2e9ee5751d56f6c6282d15d9b006f5b18eed9d))
* require nonexistent submission download directory ([#102](https://github.com/BfArM-MVH/grz-cli/issues/102)) ([a42b731](https://github.com/BfArM-MVH/grz-cli/commit/a42b7315f0fd78b3f337bfd4ace443c43b577a51))
* skip strict read length check for ONT reads ([#107](https://github.com/BfArM-MVH/grz-cli/issues/107)) ([1bd2055](https://github.com/BfArM-MVH/grz-cli/commit/1bd2055fb433b9ddb91b3ddcd77dabaabe81d777))

## [0.5.0](https://github.com/BfArM-MVH/grz-cli/compare/v0.4.0...v0.5.0) (2025-04-28)


### Features

* add 'list' command to list inbox submissions ([dbc4f6d](https://github.com/BfArM-MVH/grz-cli/commit/dbc4f6d300c6604e822bdfacc0e75cc8409d925d))
* block fully-uploaded submissions from being reuploaded ([#83](https://github.com/BfArM-MVH/grz-cli/issues/83)) ([bc82a30](https://github.com/BfArM-MVH/grz-cli/commit/bc82a303509f11e270c8d1e729047f76d0f2e29a))
* fallback to environment for S3 credentials ([#88](https://github.com/BfArM-MVH/grz-cli/issues/88)) ([637d95a](https://github.com/BfArM-MVH/grz-cli/commit/637d95acc1dcd2d970684d0642fdf4c8ce9dfff8))
* introduce `submit` subcommand for sequentially running validate, encrypt and upload ([#69](https://github.com/BfArM-MVH/grz-cli/issues/69)) ([6235a4c](https://github.com/BfArM-MVH/grz-cli/commit/6235a4cae5a485c01f2eb28f67289d1874aa8183))
* introduce grz-mode for internal use ([#74](https://github.com/BfArM-MVH/grz-cli/issues/74)) ([1b836d7](https://github.com/BfArM-MVH/grz-cli/commit/1b836d700259663e603d36a79fe48bb225a7dca5))


### Performance Improvements

* Increase minimum chunksize ([#72](https://github.com/BfArM-MVH/grz-cli/issues/72)) ([10dfec9](https://github.com/BfArM-MVH/grz-cli/commit/10dfec9369b907e2eef1a603c07793491c9e64f1))


### Documentation

* add new S3 option to example config ([#71](https://github.com/BfArM-MVH/grz-cli/issues/71)) ([4eeff53](https://github.com/BfArM-MVH/grz-cli/commit/4eeff5332a4d9a8b44886b00a2795eb403d591f4))

## [0.4.0](https://github.com/BfArM-MVH/grz-cli/compare/v0.3.0...v0.4.0) (2025-04-16)


### Features

* add BAM file validation (enables long read submissions)
* fail fast during FASTQ validation
* S3 support for special backends
* generate a submission ID as upload target folder

## [0.3.0](https://github.com/BfArM-MVH/grz-cli/compare/v0.2.0...v0.3.0) (2025-03-21)


### Features

* report supported metadata schema versions with cli version ([#49](https://github.com/BfArM-MVH/grz-cli/issues/49)) ([07d861e](https://github.com/BfArM-MVH/grz-cli/commit/07d861e32302632ca55a7f5cb5ebbfb05d4a2649))


### Bug Fixes

* use metadata model from grz-pydantic-models ([#41](https://github.com/BfArM-MVH/grz-cli/issues/41)) ([5f2944e](https://github.com/BfArM-MVH/grz-cli/commit/5f2944e7892820a9ba4e629633f8a50df4cb2a01))

## [0.2.0](https://github.com/BfArM-MVH/grz-cli/compare/v0.1.4...v0.2.0) (2025-03-18)


### Features

* Update metadata model to v1.1.1 ([#35](https://github.com/BfArM-MVH/grz-cli/issues/35)) ([eda8899](https://github.com/BfArM-MVH/grz-cli/commit/eda88993e13c62240fde8f4ea56dc7ea7c5a96e6))


### Bug Fixes

* `grz-cli --version` should report correct version ([#37](https://github.com/BfArM-MVH/grz-cli/issues/37)) ([a185e42](https://github.com/BfArM-MVH/grz-cli/commit/a185e42fc8f577598db2cdd4a007e816b6856ea8))
* decryption log message should read "decryption", not "encryption" ([#36](https://github.com/BfArM-MVH/grz-cli/issues/36)) ([70708dd](https://github.com/BfArM-MVH/grz-cli/commit/70708dd65bd342bedb27e6f5cb3f521b276bcf01))
* Ensure exact R1/R2 mappings in paired end sequencing submissions ([#22](https://github.com/BfArM-MVH/grz-cli/issues/22)) ([cdc3f9a](https://github.com/BfArM-MVH/grz-cli/commit/cdc3f9a1968950a45c20a35e39d4be504635c9c4))
* ensure reference genomes are consistent ([#26](https://github.com/BfArM-MVH/grz-cli/issues/26)) ([849c6b2](https://github.com/BfArM-MVH/grz-cli/commit/849c6b25d8a92b0d68be45c7ba31348b1f937681))

## [0.1.4](https://github.com/BfArM-MVH/grz-cli/compare/v0.1.3...v0.1.4) (2025-02-21)


### Bug Fixes

* Add check for duplicate lab data names ([#8](https://github.com/BfArM-MVH/grz-cli/issues/8)) ([ca6fc5d](https://github.com/BfArM-MVH/grz-cli/commit/ca6fc5d3ce679f3e063ec4aa2050703600a94d0c))
* Check if oncology samples have tumor cell counts. ([#7](https://github.com/BfArM-MVH/grz-cli/issues/7)) ([9341d9f](https://github.com/BfArM-MVH/grz-cli/commit/9341d9f30a3ad114d881452710f08c84fa9789ab))
* Update metadata schema ([#9](https://github.com/BfArM-MVH/grz-cli/issues/9)) ([7ca7ab0](https://github.com/BfArM-MVH/grz-cli/commit/7ca7ab0a31289706bba9217a05ba75e401d34840))

## 0.1.3 (2024-12-16)


### Bug Fixes

* check all files' existence before uploading any file at all ([47b2b91](https://github.com/BfArM-MVH/GRZ_CLI/commit/47b2b9176a7552b409c33f2f95150130d9ff3fa1))
* Custom help text for top-level grz-cli ([675e4af](https://github.com/BfArM-MVH/GRZ_CLI/commit/675e4af97a93e4b19d47ac3f13aa0f66d2ad811e))
* force endpoint_url to be AnyHttpUrl ([0ca9111](https://github.com/BfArM-MVH/GRZ_CLI/commit/0ca9111c3becf670c3aecde7afa2716f5222b4f1))


### Miscellaneous Chores

* Release-As: 0.1.3 ([a4d0c45](https://github.com/BfArM-MVH/GRZ_CLI/commit/a4d0c45e5d361338bd85da8d67d6d002e307a397))
