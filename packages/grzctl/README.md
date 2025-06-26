# grzctl

Command-line tool for internal GRZ operations.

## Running a development version

1. Install [`uv`](https://docs.astral.sh/uv)
  - An easy way is to create a Conda environment containing `uv`.
2. Clone the `grz-tools` repository locally
3. From the repository root, use `uv run grzctl <grzctl options here>`
  - Alternatively, you can use `uv run --project path/to/repo grzctl <grzctl options here>` to run it from any directory.
    This is useful if your config uses relative paths and `grzctl` must therefore be run from a specific directory.

## Running unit tests

First, ensure `uv` is installed (see above).

To run the grz-tools integration tests, run the following from the repository root:

```
uv run tox -e 3.12
```

Some packages have their own unit tests.
Run the same command above while inside a specific package directory to run that package's unit tests, if it has any.
