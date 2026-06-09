#!/usr/bin/env bash
set -euo pipefail

if [ -n "${snakemake_params[grz_private_key_passphrase]:-}" ]; then
	export C4GH_PASSPHRASE="${snakemake_params[grz_private_key_passphrase]}"
fi

db_config="${snakemake_input[db_config_path]}"
log_stdout="${snakemake_log[stdout]}"
log_stderr="${snakemake_log[stderr]}"
inbox_config_path="${snakemake_input[inbox_config_path]}"

metadata_file_path="${snakemake_input[metadata]}"
metadata_dir="$(dirname "${metadata_file_path}")"
encrypted_files_dir="${snakemake_input[encrypted_files_dir]}"
output_files_dir="${snakemake_output[files_dir]}"
mkdir -p "${output_files_dir}"
progress_logs_dir="$(dirname "${snakemake_output[progress_log]}")"

# grzctl decrypt handles DB state transitions (DECRYPTING → DECRYPTED) via DbContext.
grzctl decrypt \
	--config-file "${inbox_config_path}" \
	--config-file "${db_config}" \
	--metadata-dir "${metadata_dir}" \
	--encrypted-files-dir "${encrypted_files_dir}" \
	--files-dir "${output_files_dir}" \
	--logs-dir "${progress_logs_dir}" \
	>"$log_stdout" 2>"$log_stderr"
