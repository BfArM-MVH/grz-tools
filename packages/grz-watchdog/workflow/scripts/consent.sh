#!/usr/bin/env bash
set -euo pipefail

log_file="${snakemake_log[stderr]}"

grzctl consent --submission-dir "${snakemake_input[data]}" >"${snakemake_output[consent_flag]}" 2>"$log_file"
