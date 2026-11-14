#!/usr/bin/env bash
set -euo pipefail

grzctl_config="${snakemake_input[grzctl_config_path]}"
log_stdout="${snakemake_log[stdout]}"
log_stderr="${snakemake_log[stderr]}"

metadata_file_path="${snakemake_input[metadata]}"
re_encrypted_files_dir="${snakemake_input[re_encrypted_files_dir]}"

read -r -a progress_logs_to_archive <<<"${snakemake_input[progress_logs_to_archive]}"
progress_logs_dir="$(dirname "${progress_logs_to_archive[0]}")"

# grzctl archive derives consent from metadata and handles DB state transitions (ARCHIVING → ARCHIVED) via DbContext.
# It redacts the metadata and the logs before it uploads them.
grzctl --config "${grzctl_config}" archive \
	--metadata-dir "$(dirname "$metadata_file_path")" \
	--logs-dir "${progress_logs_dir}" \
	--encrypted-files-dir "${re_encrypted_files_dir}" \
	>"$log_stdout" 2>"$log_stderr"
