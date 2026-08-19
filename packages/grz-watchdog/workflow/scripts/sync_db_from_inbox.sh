#!/usr/bin/env bash
set -euo pipefail

# Scan the inbox and register any new/updated submissions in the database.
grzctl --config "${snakemake_input[grzctl_config_path]}" db sync-from-inbox \
	>"${snakemake_log[stdout]}" 2>"${snakemake_log[stderr]}"
