#!/usr/bin/env bash
set -euo pipefail

submitter_id="${snakemake_wildcards[submitter_id]}"
inbox="${snakemake_wildcards[inbox]}"

# Scan the inbox and register any new/updated submissions in the database.
grzctl --config "${snakemake_input[grzctl_config_path]}" db sync-from-inbox \
	--submitter-id "${submitter_id}" \
	--inbox "${inbox}" \
	>"${snakemake_log[stdout]}" 2>"${snakemake_log[stderr]}"
