#!/usr/bin/env bash
set -euo pipefail

if [ -n "${snakemake_params[s3_access_key]:-}" ]; then
	export GRZ_S3__ACCESS_KEY="${snakemake_params[s3_access_key]}"
fi

if [ -n "${snakemake_params[s3_secret]:-}" ]; then
	export GRZ_S3__SECRET="${snakemake_params[s3_secret]}"
fi

submission_id="${snakemake_wildcards[submission_id]}"
submitter_id="${snakemake_wildcards[submitter_id]}"
inbox="${snakemake_wildcards[inbox]}"
grzctl_config="${snakemake_input[grzctl_config_path]}"
log_stdout="${snakemake_log[stdout]}"
log_stderr="${snakemake_log[stderr]}"
output_metadata_dir="${snakemake_output[metadata_dir]}"
mkdir -p "${output_metadata_dir}"
output_encrypted_files_dir="${snakemake_output[encrypted_files_dir]}"
progress_logs_dir="$(dirname "${snakemake_output[progress_log]}")"
mkdir -p "${progress_logs_dir}"

# grzctl download handles DB state transitions (DOWNLOADING → DOWNLOADED) and
# submission metadata population via DbContext (--update-db is the default).
grzctl --config "${grzctl_config}" download \
	--inbox "${inbox}" \
	--submission-id "${submission_id}" \
	--metadata-dir "${output_metadata_dir}" \
	--encrypted-files-dir "${output_encrypted_files_dir}" \
	--logs-dir "${progress_logs_dir}" \
	>"$log_stdout" 2>"$log_stderr"
