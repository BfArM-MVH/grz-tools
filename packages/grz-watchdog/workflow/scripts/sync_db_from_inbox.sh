#!/usr/bin/env bash
set -euo pipefail

# Scan the inbox and register any new/updated submissions in the database.
# grzctl db sync-from-inbox:
#   - merges the inbox (S3) config for query_submissions
#   - merges the DB config (via the db group's --config-file) for DB access
#   - adds new submissions and transitions uploading → uploaded as needed
grzctl db \
	--config-file "${snakemake_input[inbox_config_path]}" \
	--config-file "${snakemake_input[db_config_path]}" \
	sync-from-inbox \
	>"${snakemake_log[stdout]}" 2>"${snakemake_log[stderr]}"

