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
Load it in your shell's startup file:

```bash
# ~/.bashrc
eval "$(grzctl completion bash)"
# ~/.zshrc
eval "$(grzctl completion zsh)"
# ~/.config/fish/config.fish
grzctl completion fish | source
```

This runs `grzctl` at every shell start.
To avoid that, save the output to a file once and source that file instead.

Set `GRZCTL_COMPLETE_SUBMISSION_IDS=1` to also complete submission IDs in `grzctl db` from the database.
Each completion then loads the config and queries the database.

## Documentation

- [Registering change requests](docs/change-requests.md) — recording that a submission must be modified, deleted, or transferred, with the who/when/what audit trail.
- [Crypt4GH keys](docs/crypt4gh-keys.md): which Crypt4GH keys grz-cli and grzctl need, where they go in each config, and the key layouts of grzctl.
