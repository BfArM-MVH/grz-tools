# `grz-cli`

A command-line tool for validating, encrypting, and uploading submissions to a genomDE Model Project GDC (Genome Data Center).


## Table of Contents

- [Installation](#installation)
  - [Requirements](#requirements)
  - [Using Conda (recommended)](#using-conda-recommended)
  - [Using Docker](#using-docker)
- [Usage](#usage)
  - [Configuration](#configuration)
  - [Submission Layout](#submission-layout)
  - [Example submission procedure](#example-submission-procedure)
  - [Troubleshooting](#troubleshooting)
- [Command-Line Interface](#command-line-interface)
- [Contributing](#contributing)
- [License](#license)
- [Acknowledgements](#acknowledgements)


## Installation

### Requirements

The current version of the tool requires the working directory to have at least as much free disk space as the total size of the data being submitted.

Beside of the disk space requirements for the submission data, this tool also requires a Linux environment.
For example:

- Linux server
- Virtual machine running Linux
- Docker container
- Windows Subsystem for Linux


### Using [Conda](https://conda.io) (recommended)

If Conda is not yet available on your system, we recommend to install it through the [Miniforge Conda installer](https://github.com/conda-forge/miniforge) by running the following commands:

```bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
```

Next, install the `grz-cli` tool:

```bash
conda create -n grz-tools -c conda-forge -c bioconda grz-cli
conda activate grz-tools
grz-cli --help
```

#### Updating

Use the following command to update the tool:

```bash
conda update -n grz-tools -c conda-forge -c bioconda grz-cli
```


### Using pip (not recommended)

While installation via `pip` is possible, it is not recommended because users should create/manage a Python virtual environment and must ensure that the correct Python version is being used.

```bash
pip install grz-cli
```

#### Updating

Use the following command to update the tool:

```bash
pip upgrade grz-cli
```


### Using Docker

Docker images are available via biocontainers at [https://biocontainers.pro/tools/grz-cli](https://biocontainers.pro/tools/grz-cli).

The build process can take at least a few days after the Bioconda release, so double-check that the [latest version in Bioconda](https://anaconda.org/bioconda/grz-cli) is also the latest Docker image version.

## Usage

### Configuration

**The configuration file will be provided by your associated GRZ.**
**Do not create this file as an LE.**

The tool requires a configuration file in YAML format to specify the S3 API parameters and other validation options.

This file may be placed at `~/.config/grz-cli/config.yaml` or provided each time to `grz-cli` using the `--config-file` option on the command line.

The S3 secrets can either be directly within the config file or as defined with the usual AWS environment variables: `AWS_ACCESS_KEY_ID` and `AWS_SECRET_ACCESS_KEY`.

### Submission Layout

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
│   └─── target_regions.bed
└── metadata
    └── metadata.json
```

The only requirements are that `metadata/metadata.json` exists and the `files/` directory contains all of the other files.
Data files may be nested under subfolders inside `files/` for better organization.
For example, each donor could have their own folder for files.


### Example submission procedure

After preparing your submission as outlined above, you can use the following command to validate, encrypt and upload the submission:

```bash
grz-cli submit --submission-dir EXAMPLE_SUBMISSION
```

### Troubleshooting

In case of issues, please re-run your commands with `grz-cli --log-level DEBUG --log-file path/to/write/file.log [...]` and submit the log file to the GDC data steward.

## Command-Line Interface

`grz-cli` provides a command-line interface with the following subcommands:

### submit

The `submit` command is the recommended command for submitting data to genome data centers.
It combines the `validate`, `encrypt`, and `upload` commands (see below).

- `-s, --submission-dir`: Path to the submission directory containing both 'metadata/' and 'files/' directories [**Required**]
- `-c, --config-file`: Path to config file [_optional_]

```bash
grz-cli submit --submission-dir foo
```


### validate

It is recommended to run this command before continuing with encryption and upload.

- `--submission-dir`: Path to the submission directory containing both 'metadata/' and 'files/' directories [**Required**]

Example usage:

```bash
grz_cli validate --submission-dir foo
```

### encrypt

If a working directory is not provided, then the current directory is used automatically.
Files are stored in a folder named `encrypted_files` as a sub-folder of the working directory.

- `-s, --submission-dir`: Path to the submission directory containing both 'metadata/' and 'files/' directories [**Required**]
- `-c, --config-file`: Path to config file [_optional_]

```bash
grz-cli encrypt --submission-dir foo
```

### upload

Upload the submission into a S3 structure of a GRZ.

- `-s, --submission-dir`: Path to the submission directory containing both `metadata/` and `encrypted_files/` directories [**Required**]
- `-c, --config-file`: Path to config file [_optional_]

Example usage:

```bash
grz-cli upload --submission-dir foo
```

### get-id

*Available in grz-cli v1.2.0 or higher.*

Compute and print the submission ID from a submission's JSON metadata.

This is useful in case you forget to store the ID printed during upload.

Example usage:

```bash
grz-cli get-id path/to/metadata.json
```

## Contributing

### Running unreleased/development versions

First, install `uv`.
We recommend using Conda or [Pixi](https://pixi.sh/latest).

After cloning the desired branch of the `grz-tools` repo locally, you can run `grz-cli` directly from the repo using:

```
uv run --project path/to/cloned/grz-tools grz-cli --help
```

### Testing

Please note that binary files used for testing are managed with [Git LFS](https://git-lfs.com), which will be needed to clone them locally with the git repository.

To run the tests, navigate to the root directory of your project and invoke `pytest`.
Alternatively, install `uv` and `tox` and run `uv run tox`.

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

## Acknowledgements

Parts of [Crypt4GH](https://github.com/EGA-archive/crypt4gh) are used in modified form.
