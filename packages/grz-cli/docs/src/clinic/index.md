# Instructions for participating clinics

## Installation

### Requirements

Linux is the only supported operating system.
A Linux environment can be acquired in many different ways, for example:

- A server running Linux directly
- A virtual machine running Linux
- A Docker container
- Windows Subsystem for Linux

Currently, the encrypted copies of each file are cached within the submission folder.
This means that at least the total size of the submission in extra free disk space is needed before starting.

### Using [Conda](https://conda.io) (recommended)

If Conda is not yet available on your system, we recommend to install it through the [Miniforge Conda installer](https://github.com/conda-forge/miniforge) by running the following commands:

```bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash "Miniforge3-$(uname)-$(uname -m).sh"
```

Next, install `grz-cli`:

```bash
conda create -n grz-tools -c conda-forge -c bioconda grz-cli
conda activate grz-tools
grz-cli --help
```

#### Updating

Use the following command to update `grz-cli`:

```bash
conda update -n grz-tools -c conda-forge -c bioconda grz-cli
```

### Using Docker

Docker images are unofficially available from [BioContainers](https://biocontainers.pro/tools/grz-cli).
These are automatically built from the Bioconda package and the `grz-cli` developers do not control their availability or tags.

The build process can take at least a few days after the Bioconda release, so please double-check that the latest Docker image version available from BioContainers is also the [latest version in Bioconda](https://anaconda.org/bioconda/grz-cli).

### Using pip

While installation via `pip` is possible, it is **not recommended** because it requires careful consideration of the desired installation environment and ensuring that a supported Python version is being used.

```bash
pip install grz-cli
```

#### Updating

Use the following command to update `grz-cli`:

```bash
pip upgrade grz-cli
```

## Usage

### Configuration

!!! info
    The configuration file will be provided by your associated GDC.
    Do not create this file yourself.

The tool requires a configuration file in YAML format to specify your genome data center's S3 API parameters, their inbox public key for encryption, and other validation options.

This file may be placed at `~/.config/grz-cli/config.yaml` or provided each time to `grz-cli` using the `--config-file` option on the command line.

The S3 secrets can either be directly defined within the config file or with the usual AWS environment variables: `AWS_ACCESS_KEY_ID` and `AWS_SECRET_ACCESS_KEY`.

### Submission Layout

It is recommended to have the following folder structure for a single submission:

```
EXAMPLE_SUBMISSION
├── files
│   ├── donor1_blood_normal.read1.fastq.gz
│   ├── donor1_blood_normal.read2.fastq.gz
│   ├── donor1_blood_normal.vcf
│   ├── donor1_blood_tumor.read1.fastq.gz
│   ├── donor1_blood_tumor.read2.fastq.gz
│   ├── donor1_blood_tumor.vcf
│   └─── target_regions.bed
└── metadata
    └── metadata.json
```

!!! warning
    Do not use the tanG anywhere in your file names.
    GDCs cannot keep the tanG long-term and file names are archived without modification.

The only requirements are that `metadata/metadata.json` exists and the `files/` directory contains all of the other files.
Data files can be nested under subfolders inside `files/` for better organization.
For example, each donor could have their own folder for files.

### Submitting

After preparing your submission as outlined above, you can use the following command to validate, encrypt, and upload the submission all at once:

```bash
grz-cli submit --submission-dir EXAMPLE_SUBMISSION/
```

Alternatively, there are separate subcommands for each step.
See `grz-cli --help` or [this page](../cli.md) for more information on the command line interface of `grz-cli`.

### Troubleshooting

In case of issues, please re-run your commands with `grz-cli --log-level DEBUG --log-file path/to/write/file.log [...]` and submit the log file to your GDC's data steward.
