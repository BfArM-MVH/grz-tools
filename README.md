# grz-tools

This is the monorepo for the project-specific software used to interact with and administer a genome data center in the genomeDE Model Project.
It contains the following packages:

- [`grz-cli`](packages/grz-cli/README.md) - A command-line tool for validating, encrypting and uploading submissions to a GRZ.
- [`grz-pydantic-models`](packages/grz-pydantic-models/README.md) - Pydantic models for schemas related to the genomDE Model Project.
- [`grzctl`](packages/grzctl/README.md) - GRZ internal tooling.
- [`grz-common`](packages/grz-common/README.md) - Common code shared between packages in `grz-tools`.
- [`grz-db`](packages/grz-db/README.md) - Libraries, SQL models and alembic migrations for the GRZ internal submission DB.

**If you are staff at a participating clinic (LE)**, you should start with the [`grz-cli` documentation](https://BfArM-MVH.github.io/grz-tools/grz-cli/).

## Contributing

GitHub Actions runs a number of quality checks for each pull request.
To minimize the number of errors triggered by your code changes, run at least the following and make sure it all passes without error before requesting a review:

```
uv run ruff format
uv run ruff check --fix
uv run tox -e typecheck
uv run tox -e 3.13
```

### Testing

Please note that binary files used for testing are managed with [Git LFS](https://git-lfs.com), which will be needed to clone them locally with the git repository.

To run the workspace-level unit tests (those in the top-level `tests/` folder) against Python 3.13, run the following:

```
uv run tox -e 3.13
```

To run the unit tests of a specific subpackage, do the same from that package's root.
For example:

```
cd packages/grz-pydantic-models
uv run tox -e 3.13
```
