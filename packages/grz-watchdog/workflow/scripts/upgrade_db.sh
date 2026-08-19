#!/usr/bin/env bash
set -euo pipefail

grzctl --config "${snakemake_input[grzctl_config_path]}" db upgrade >"${snakemake_log[stdout]}" 2>"${snakemake_log[stderr]}"
