#!/bin/bash
set -euo pipefail

pixi run --manifest-path /grz-watchdog/pixi.toml -- \
  snakemake \
    --workflow-profile /grz-watchdog/workflow/profiles/default \
    --snakefile /grz-watchdog/workflow/Snakefile \
    --directory /workdir \
    --conda-prefix /conda-envs \
    "$@"
