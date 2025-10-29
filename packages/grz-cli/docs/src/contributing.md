# Contributing

## Running unreleased/development versions

First, install `uv`.
We recommend using Conda or [Pixi](https://pixi.sh/latest).

After cloning the desired branch of the `grz-tools` repo locally, you can run `grz-cli` directly from the cloned source repository using:

```
uv run --project path/to/cloned/grz-tools grz-cli --help
```

## Building the documentation

```
cd packages/grz-cli
uv run mike deploy --config-file docs/mkdocs.yaml development
uv run mike serve --config-file docs/mkdocs.yaml
```

This will build the `grz-cli` docs, deploy them to your local `gh-pages` branch, and start a webserver serving from that local branch.

Use the version selection drop-down in the header bar to load the "development" docs you just deployed instead of the default stable docs.
