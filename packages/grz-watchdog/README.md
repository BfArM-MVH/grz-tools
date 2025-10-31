# grz-watchdog

A prototype of a snakemake workflow for the grz-ingest process.
Can be run to process a single (known) target (see `process-single`), a batch of (unknown) targets (see `process-batch`)
or continuously scan inboxes for new targets (see `daemon`).

## Setup

### Installation

This project is managed using `pixi` and requires a recent pixi version (`>=0.58.0`).

```sh
pixi install
```

For faster validation, also install `grz-check` into the same environment (not installed by default).

Alternatively, build a docker image from the provided `Dockerfile` and run that – this does include grz-check by
default. Its default entrypoint is `daemon`.

### Configuration

The following environment variables allow customization of certain paths:

```sh
GRZ_WATCHDOG_WORKDIR=/path/to/workdir  # where watchdog will store its results and basepath for relative dirs.
GRZ_WATCHDOG_SNAKEMAKE_PROFILE=/path/to/your/custom/snakemake/profile
```

Setting `GRZ_WATCHDOG_WORKDIR` is strongly recommended, otherwise the default workdir is relative to the project
directory of grz-watchdog, i.e., `grz-watchdog/workdir`.

Example grz-watchdog configuration (refer to schema at `workflow/schemas/config.schema.yaml` for details):

```yaml
monitor:
  interval: 30  # The interval in seconds for monitoring.
batch:
  limit: 5  # The batch processing limit.

estimates: # Estimates of certain properties of grz-watchdog and its components
  speed: # Estimates of download or processing speeds in MB/s
    download: 200
    decrypt: 250
    validate: 100
    encrypt: 250
    archive: 200
    qc: 10

temp-outputs: true  # If true, rule outputs are temporary (marked as `temp(…)`). If false, rule outputs will be kept.
on-failed-validation: "cleanup"  # Action to take on failed validation.
auto-cleanup: "inbox+storage"  # "The auto-cleanup strategy."

handlers:
  # What to do when the workflow encounters an error.
  # Placeholder '{log}' points to main snakemake log.
  # For example, to send a mail with the contents of the log:
  #   "mail -s 'grz-watchdog error' send-mail-to@recipient.domain < {log}"
  on-error: "echo 'Failed running grz-watchdog, check logs at {log}'; cat {log}"
  on-start: "echo 'Starting grz-watchdog.'"
  on-success: "echo 'Successfully finished grz-watchdog run.'"


config_paths: # Paths to grzctl config files. Prefer absolute paths.
  inbox:
    "123456789": # le id
      inbox1: "/path/to/configs/inbox.le1.inbox1.yaml"
      inbox2: "/path/to/configs/inbox.le1.inbox2.yaml"  # if you have more than 1 inbox per le
    "234567890": # yet another le id
      inbox1: "/path/to/configs/inbox.le2.inbox1.yaml"
      # … etc …
  archive:
    consented: "/path/to/configs/consented.yaml"
    nonconsented: "/path/to/configs/nonconsented.yaml"
  db: "/path/to/configs/db.yaml"
  pruefbericht: "/path/to/configs/pruefbericht.yaml"

qc:
  revision: "v1.2.0"
  reference_directory: "/path/to/qc/references"
  selection_strategy:
    enabled: false  # Whether automatic qc selection is enabled or not.
    target_percentage: 2.0  # Target percentage of submissions to run QC on per LE per month.
    salt: 鹽  # Salt for random seed
  run-qc:
    extra: "" # extra args to GRZ_QC_Workflow nextflow invocation
```

### QC references

If you already have reference indices for the QC workflow, copy them to the `reference_directory` specified in the
config and `mkdir -p "<resources>/shared_qc_launchdir"` (a marker so snakemake doesn't unnecessarily re-run building
reference indices).
Otherwise, `grz-watchdog` will build the references the first time they are needed.

### Resources

To help snakemake schedule jobs, define available resources either in your snakemake profile or 
via the corresponding cli options.
For example
```yaml
cores: 64
resources:
  mem_mb: 60_000      # total available memory
  disk_mb: 4_000_000  # total available storage space
```

### Miscellaneous remarks

- to control where pixi environments are stored, refer
  to [detached-environments](https://pixi.sh/dev/reference/pixi_configuration/#detached-environments)
- to control where pixi packages caches are stored, refer to TODO
- to control where nextflow's conda caches are stored, refer to TODO

## Examples

A general remark: Any non grz-watchdog specific arguments/options are passed as is on to snakemake.
That means you can easily:

- specify snakemake profiles to use,
- enable verbose output (`-v`/`-vv`/`…`),
- use custom loggers,
- enable `--keep-going`
- enable `--rerun-incomplete`,
- control how temp files are handled,
- … etc …

### Processing a specific submission ID

Process a specific submission by supplying the submitter ID, inbox name and submission ID.
By default, whether QC will be performed is determined automatically using the selection strategy defined in the
grz-watchdog config:

```sh
SUBMITTER_ID="123456789"
INBOX="inbox1"
SUBMISSION_ID="123456789_2025-09-22_aabbccdd"
pixi run process-single ${SUBMITTER_ID} ${INBOX} ${SUBMISSION_ID} --configfile /path/to/watchdog/config.yaml
```

To forcefully enable/disable QC, use options `--force-with-qc/--force-without-qc`, e.g.:

```sh
pixi run process-single ${SUBMITTER_ID} ${INBOX} ${SUBMISSION_ID} --force-without-qc --configfile /path/to/watchdog/config.yaml
```

### Batch Processing

Depending on the `batch: limit` in `config/config.yaml`, determine and take up to `limit` submissions ready for
processing (by means of an automatic inbox scan) and submit these:

```sh
pixi run process-batch --configfile /path/to/watchdog/config.yaml
```

(will be executed within `$GRZ_WATCHDOG_WORKDIR` by default; can be overwritten by appending
`--directory /some/other/path` or re-defining `$$GRZ_WATCHDOG_WORKDIR`).

### Daemon Mode

To start the daemon, invoke:

```sh
pixi run daemon --configfile /path/to/watchdog/config.yaml
```

This will constantly monitor the configured inboxes and check for new submissions and sync them to the database.
Then, submissions available for processing (determined by their state in the database) are pushed to the job queue
automatically.