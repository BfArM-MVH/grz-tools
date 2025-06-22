# GRZ Tools Monorepo

This monorepo hosts the following packages:

- [`grz-cli`](packages/grz-cli/README.md) - A command-line tool for validating, encrypting and uploading submissions to a GRZ.
- [`grz-pydantic-models`](packages/grz-pydantic-models/README.md) - Pydantic models for schemas related to the genomDE Model Project.
- [`grzctl`](packages/grzctl/README.md) - GRZ internal tooling.
- [`grz-common`](packages/grz-common/README.md) - Common code shared between packages in `grz-tools`.
- [`grz-db`](packages/grz-db/README.md) - Libraries, SQL models and alembic migrations for the GRZ internal submission DB.

## grz-cli

The [`grz-cli`](packages/grz-cli/README.md) package is the primary CLI for submissions to the GRZs.
It provides functionality for:
- Validating submissions
- Encrypting files using crypt4gh
- Uploading files to a GRZ

For detailed installation and usage instructions, please refer to the [grz-cli README](packages/grz-cli/README.md).

## Legacy Information

Previous grz-cli repository content is still available in the `archive/pre-monorepo` branch.
