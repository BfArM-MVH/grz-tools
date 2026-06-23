#!/usr/bin/env bash
set -euo pipefail

if [ -n "${snakemake_params[s3_access_key]:-}" ]; then
	export GRZ_S3__ACCESS_KEY="${snakemake_params[s3_access_key]}"
fi

if [ -n "${snakemake_params[s3_secret]:-}" ]; then
	export GRZ_S3__SECRET="${snakemake_params[s3_secret]}"
fi

submission_id="${snakemake_wildcards[submission_id]}"
db_config="${snakemake_input[db_config_path]}"
log_stdout="${snakemake_log[stdout]}"
log_stderr="${snakemake_log[stderr]}"
output_metadata_dir="${snakemake_output[metadata_dir]}"
mkdir -p "${output_metadata_dir}"
output_encrypted_files_dir="${snakemake_output[encrypted_files_dir]}"
progress_logs_dir="$(dirname "${snakemake_output[progress_log]}")"
mkdir -p "${progress_logs_dir}"

# grzctl download handles DB state transitions (DOWNLOADING → DOWNLOADED) and
# submission metadata population via DbContext (--update-db is the default).
grzctl download \
	--config-file "${snakemake_input[inbox_config_path]}" \
	--config-file "${db_config}" \
	--submission-id "${submission_id}" \
	--metadata-dir "${output_metadata_dir}" \
	--encrypted-files-dir "${output_encrypted_files_dir}" \
	--logs-dir "${progress_logs_dir}" \
	>"$log_stdout" 2>"$log_stderr"
