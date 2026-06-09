#!/usr/bin/env bash
set -euo pipefail

db_config="${snakemake_input[db_config_path]}"
log_stdout="${snakemake_log[stdout]}"
log_stderr="${snakemake_log[stderr]}"

metadata_file_path="${snakemake_input[metadata]}"
metadata_dir="$(dirname "$metadata_file_path")"
unencrypted_files_dir="${snakemake_input[files_dir]}"
output_encrypted_files_dir="${snakemake_output[encrypted_files_dir]}"
mkdir -p "${output_encrypted_files_dir}"
progress_logs_dir="$(dirname "${snakemake_output[encryption_log]}")"

CONSENT=$(cat "${snakemake_input[consent_flag]}")
if [[ "$CONSENT" == "true" ]]; then
	CONFIG_FILE="${snakemake_input[consented_config_path]}"
else
	CONFIG_FILE="${snakemake_input[nonconsented_config_path]}"
fi

echo "Consent: $CONSENT. Using config file: $CONFIG_FILE" >"$log_stdout" 2>"$log_stderr"

# grzctl encrypt handles DB state transitions (ENCRYPTING → ENCRYPTED) via DbContext.
grzctl encrypt \
	--config-file "$CONFIG_FILE" \
	--config-file "${db_config}" \
	--metadata-dir "${metadata_dir}" \
	--files-dir "${unencrypted_files_dir}" \
	--output-encrypted-files-dir "${output_encrypted_files_dir}" \
	--logs-dir "${progress_logs_dir}" \
	--force \
	>>"$log_stdout" 2>>"$log_stderr"

