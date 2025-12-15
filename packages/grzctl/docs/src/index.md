!!! warning "Only relevant for GDC staff"
    If you are staff at a clinic looking to upload submissions to a genome data center, please see the [`grz-cli` docs](https://BfArM-MVH.github.io/grz-tools/grz-cli/) instead.

Welcome to the documentation for `grzctl`, a command-line tool for managing internal GRZ operations.

## Running a development version

1. Install [`uv`](https://docs.astral.sh/uv).
    - An easy way is to create a Conda environment containing `uv`.
2. Clone the `grz-tools` repository locally.
3. From the repository root, use `uv run grzctl <grzctl options here>`.
    - Alternatively, you can use `uv run --project path/to/repo grzctl <grzctl options here>` to run it from any directory.
    This is useful if your config uses relative paths and `grzctl` must therefore be run from a specific directory.

## Running unit tests

First, ensure `uv` is installed (see above).

``` shell
cd packages/grzctl
uv run tox -e 3.12
```

## Debugging Textual

Start the remote Textual debugging console on the machine you will use to debug the Textual app.

``` shell
uv run textual console
```

Use the `TEXTUAL=devtools` environment variable to instruct Textual to connect to the remote debug console.

``` shell
TEXTUAL=devtools uv run grzctl db --config-file config.db.yaml tui
```
