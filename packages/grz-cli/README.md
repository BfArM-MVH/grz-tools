# GRZ CLI

A command-line tool for validating, encrypting, uploading and downloading submissions to/from a GDC/GRZ (Genomrechenzentrum).

## Table of Contents

- [Introduction](#introduction)
- [Features](#features)
- [Installation](#installation)
  - [Requirements](#requirements)
  - [End-user setup](#end-user-setup)
    - [Installation via `conda` (recommended)](#installation-via-conda-recommended)
    - [Installation via `pip` (not recommended)](#installation-via-pip-not-recommended)
- [Usage](#usage)
  - [Configuration](#configuration)
  - [Exemplary submission procedure](#exemplary-submission-procedure)
- [Command-Line Interface](#command-line-interface)
- [Testing](#testing)
- [Contributing](#contributing)
- [License](#license)
- [Acknowledgements](#acknowledgements)

## Introduction

This tool provides a way to validate files, encrypt/decrypt files using the [crypt4gh](https://crypt4gh.readthedocs.io/en/latest/) library and upload/download the encrypted files to an S3 bucket of a GDC/GRZ. It also logs the progress and outcomes of these operations in a metadata file.

It is recommended to have the following folder structure for a single submission:

```
EXAMPLE_SUBMISSION
├── files
│   ├── aaaaaaaa00000000aaaaaaaa00000000_blood_normal.read1.fastq.gz
│   ├── aaaaaaaa00000000aaaaaaaa00000000_blood_normal.read2.fastq.gz
│   ├── aaaaaaaa00000000aaaaaaaa00000000_blood_normal.vcf
│   ├── aaaaaaaa00000000aaaaaaaa00000000_blood_tumor.read1.fastq.gz
│   ├── aaaaaaaa00000000aaaaaaaa00000000_blood_tumor.read2.fastq.gz
│   ├── aaaaaaaa00000000aaaaaaaa00000000_blood_tumor.vcf
│   ├── target_regions.bed
└── metadata
    └── metadata.json
```

The current version of the tool requires the `working_dir` to have at least as much free disk space as the total size of the data being submitted.
## Features

- **Validation**: Validate file checksums, basic file metadata and BfArM requirements.
- **Encryption**: Encrypt files using `crypt4gh`.
- **Decryption**: Encrypt files using `crypt4gh`.
- **Upload**: Upload encrypted files directly to a GRZ either (via built-in `boto3`).
- **Download**: Download encrypted files from a GRZ (via built-in `boto3`).
- **Logging**: Log progress and results of operations


## Installation

### Requirements
Beside of the disk space requirements for the submission data, this tool also requires a linux environment, e.g.:
- Linux server
- Virtual machine running linux
- Docker container
- Windows subsystem for linux
- ...


### End-user setup
The recommended method to install this tool is using the conda package manager.


#### Installation via `conda` (recommended)

If `conda` is not yet available on your system, we recommend to install the [Miniforge conda distribution](https://github.com/conda-forge/miniforge) by running the following commands:
```bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
```
There are also alternative ways to install conda:
- [Micromamba, a single executable that does not require a base environment](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html)
- [Official installation instructions](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)

Next, install the `grz-cli` tool:
```bash
# create conda environment and activate it
conda create -n grz-tools -c conda-forge -c bioconda "grz-cli"
conda activate grz-tools
```

##### Update instructions

Use the following command to update the tool:

```bash
conda update -n grz-tools -c conda-forge -c bioconda grz-cli
```


#### Installation via `pip` (not recommended)
While installation via `pip` is possible, it is not recommended because users must ensure
that the correct Python version is already installed and that they are using a virtual python environment.

```bash
pip install grz-cli
```
##### Update instructions:
Use the following command to update the tool:
```bash
pip upgrade grz-cli
```


#### Docker
Docker images are available via biocontainers at [https://biocontainers.pro/tools/grz-cli](https://biocontainers.pro/tools/grz-cli).

The build process can take at least a few days after the Bioconda release, so double-check that the latest version in Bioconda is also the latest Docker image version.

## Usage

### Configuration
**The configuration file will be provided by your associated GRZ, please place it into `~/.config/grz-cli/config.yaml`.**

The tool requires a configuration file in YAML format to specify the S3 bucket and other options.
For an exemplary configuration, see [resources/config.yaml](resources/config.yaml).

S3 access and secret key can be listed either in the config file or as environment variable (`AWS_ACCESS_KEY_ID`, `AWS_SECRET_ACCESS_KEY`).

### Exemplary submission procedure
After preparing your submission as outlined above, you can use the following commands to validate, encrypt and upload the submission: 
```sh
# Validate the submission
grz-cli validate --submission-dir EXAMPLE_SUBMISSION

# Encrypt the submission
grz-cli encrypt --submission-dir EXAMPLE_SUBMISSION

# Upload the submission
grz-cli upload --submission-dir EXAMPLE_SUBMISSION
```

### Troubleshooting
**In case of issues, please re-run your commands with `grz-cli --log-level DEBUG --log-file <your-log-file.log> [...]` and submit the log file to the GRZ data steward!**

## Command-Line Interface

`grz-cli` provides a command-line interface with the following subcommands:

### validate

It is recommended to run this command before continuing with encryption and upload.
Progress files are stored relative to the submission directory.

- `--submission-dir`: Path to the submission directory containing both 'metadata/' and 'files/' directories [**Required**]

Example usage:

```bash
grz_cli validate --submission-dir foo
```

### encrypt

If a working directory is not provided, then the current directory is used automatically. The log-files are going to be stored in the sub-folder of the working directory.
Files are stored in a folder named `encrypted_files` as a sub-folder of the working directory.

- `-s, --submission-dir`: Path to the submission directory containing both 'metadata/' and 'files/' directories [**Required**]
- `-c, --config-file`: Path to config file [_optional_]

```bash
grz-cli encrypt --submission-dir foo
```

### upload

Upload the submission into a S3 structure of a GRZ.

- `-s, --submission-dir`: Path to the submission directory containing both 'metadata/' and 'encrypted_files/' directories [**Required**]
- `-c, --config-file`: Path to config file [_optional_]

Example usage:

```bash
grz-cli upload --submission-dir foo
```

## Testing

Please note that binary files used for testing are managed with [Git LFS](https://git-lfs.com), which will be needed to clone them locally with the git repository.

To run the tests, navigate to the root directory of your project and invoke `pytest`.
Alternatively, install `uv` and `tox` and run `uv run tox`.

## Contributing

<!-- Add details about how others can contribute to the project -->

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

## Acknowledgements

Parts of `crypt4gh` code is used in modified form

