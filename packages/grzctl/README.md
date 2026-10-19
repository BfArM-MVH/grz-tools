# grzctl

Command-line tool for internal GRZ operations.

## Running a development version

1. Install [`uv`](https://docs.astral.sh/uv)
  - An easy way is to create a Conda environment containing `uv`.
2. Clone the `grz-tools` repository locally
3. From the repository root, use `uv run grzctl <grzctl options here>`
  - Alternatively, you can use `uv run --project path/to/repo grzctl <grzctl options here>` to run it from any directory.
    This is useful if your config uses relative paths and `grzctl` must therefore be run from a specific directory.

## Shell completion

`grzctl completion <shell>` prints the completion script for bash (4.4 or newer), zsh or fish.
Save it once to where your shell loads completions from:

```bash
grzctl completion bash > ~/.local/share/bash-completion/completions/grzctl
grzctl completion zsh > <dir on $fpath>/_grzctl
grzctl completion fish > ~/.config/fish/completions/grzctl.fish
```

Completion needs `grzctl` on your `PATH`, so it does not work for `uv run grzctl`.
For a development version, start a shell with `uv run --project path/to/repo $SHELL`, or run `source .venv/bin/activate`.

Set `GRZCTL_COMPLETE_SUBMISSION_IDS=1` to also complete submission IDs in `grzctl db` from the database.
Each completion then loads the config and queries the database.

## Documentation

- [Registering change requests](docs/change-requests.md) — recording that a submission must be modified, deleted, or transferred, with the who/when/what audit trail.
- [Crypt4GH keys](docs/crypt4gh-keys.md): which Crypt4GH keys grz-cli and grzctl need, where they go in each config, and the key layouts of grzctl.
