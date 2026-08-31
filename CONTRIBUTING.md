# Contributing to grz-tools

This document provides guidelines for contributing to the `grz-tools` monorepo.
Please follow below steps to set up your development environment and contribute effectively.

Please also refer to the individual package READMEs for any package-specific development setup instructions.

## Development setup

- This monorepo uses [uv](https://docs.astral.sh/uv/) to manage Python virtual environments for development.
- The `grz-check` package additionally requires the [Rust toolchain](https://www.rust-lang.org/tools/install) to be
  installed.
- The `grz-db` package additionally requires the [PostgreSQL](https://www.postgresql.org/download/) database server to
  be installed.

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
uv run tox
```

Some packages have their own unit tests.
Run `uv run tox` while inside a specific package directory to run that package's unit tests.

### Database tests

Database tests should use the fixtures from `grz_db.testing` (`db` is the usual entry point) rather
than constructing their own engine. These fixtures parametrize each test over both sqlite and
PostgreSQL, so a test that builds its own engine instead only ever exercises sqlite.

The PostgreSQL half is skipped when `pg_config` is not on your `PATH`, and pytest reports this as
skips rather than failures. A green run with a large skip count can therefore mean half the database
tests never executed, so check the skip count rather than only the exit status.

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

## Docstrings

Docstrings use ReST field lists. Do not use Google-style `Args:` / `Returns:` sections.

```python
def get_submission(self, submission_id: str) -> Submission | None:
    """Retrieve a submission and its state history.

    :param submission_id: Submission ID of the submission to retrieve.
    :returns: The :class:`Submission`, or ``None`` if no submission has that ID.
    """
```

Keep the whole description above the field list. A paragraph placed after `:raises:` is rendered
detached from the description.

Mark code with double backticks. A single backtick is a different role in ReST, so `` `like this` ``
does not render as a literal.

## Static type checking

This project uses mypy for static type checking.
To run type checking, run the following from the repository root:

```bash
uv run mypy
```

## Updating dependencies

Dependencies are specified in each package's `pyproject.toml` file.
After changing dependencies, use `uv sync` to update the virtual environment for a specific package.

`uv` will only update the environment if the dependencies cannot be satisfied with the existing packages in the
environment (see also `uv.lock` file).
To update dependencies to their latest versions, use the `--upgrade` flag:

```bash
uv sync --all-packages --all-groups --all-extras --upgrade
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

