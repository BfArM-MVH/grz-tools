# Contributing to grz-tools

This document provides guidelines for contributing to the `grz-tools` monorepo.
Please follow below steps to set up your development environment and contribute effectively.

Please also refer to the individual package READMEs for any package-specific development setup instructions.

## Development setup

- This monorepo uses [uv](https://docs.astral.sh/uv/) to manage Python virtual environments for development.
- The `grz-check` package additionally requires the [Rust toolchain](https://www.rust-lang.org/tools/install) to be installed.
- The `grz-db` package additionally requires the [PostgreSQL](https://www.postgresql.org/download/) database server to be installed.

The simplest way to set up a development environment is to setup a conda environment with these dependencies:
You can use either `conda`, `mamba` or `micromamba` for this, e.g.:
```bash
mamba env create -n grz-tools -f environment-dev.yaml
mamba activate grz-tools
```

Next, install the virtual environment using `uv`:
```bash
uv sync --all-packages --all-groups --all-extras
```

## Running integration tests

To run integration tests for all packages in this monorepo, run the following from the repository root:

```bash
for d in packages/*; do
   if [ -f "$d/pyproject.toml" ]; then
      printf '\n\033[1;97;44m  RUNNING TESTS: %s  \033[0m\n\n' "$d"
      (cd "$d" && uv run tox)
   fi
done

printf '\n\033[1;97;44m  RUNNING GLOBAL TESTS  \033[0m\n\n'
uv run tox
```

Some packages have their own unit tests.
Run `uv run tox` while inside a specific package directory to run that package's unit tests.

## Code formatting and linting
This project uses ruff for code formatting and linting.

To check code formatting and linting, run the following from the repository root:

```bash
uv run ruff check
```

Some errors can automatically be fixed by ruff:

```bash
uv run ruff check --fix
```

To auto-format the code, run:

```bash
uv run ruff format
```

## Static type checking
This project uses mypy for static type checking.
To run type checking, run the following from the repository root:

```bash
uv run mypy
```

## Debugging Textual

Start the remote Textual debugging console on the machine you will use to debug the Textual app.

```
uv run textual console
```

Use the `TEXTUAL=devtools` environment variable to instruct Textual to connect to the remote debug console.

``` shell
TEXTUAL=devtools uv run grzctl db --config-file config.db.yaml tui
```

