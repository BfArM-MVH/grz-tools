#!/usr/bin/env bash
set -euo pipefail

grzctl db --config-file "${snakemake_input[db_config_path]}" upgrade >"${snakemake_log[stdout]}" 2>"${snakemake_log[stderr]}"
