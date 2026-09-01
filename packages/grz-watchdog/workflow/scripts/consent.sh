#!/usr/bin/env bash
set -euo pipefail

log_file="${snakemake_log[stderr]}"

grzctl consent --metadata-file "${snakemake_input[metadata]}" >"${snakemake_output[consent_flag]}" 2>"$log_file"
