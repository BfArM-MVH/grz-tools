#!/usr/bin/env bash
set -euo pipefail

grzctl_config="${snakemake_input[grzctl_config_path]}"
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
	ARCHIVE_FLAG="--consented"
else
	ARCHIVE_FLAG="--non-consented"
fi

echo "Consent: $CONSENT. Using archive flag: $ARCHIVE_FLAG" >"$log_stdout" 2>"$log_stderr"

# grzctl encrypt handles DB state transitions (ENCRYPTING → ENCRYPTED) via DbContext.
grzctl --config "${grzctl_config}" encrypt \
	${ARCHIVE_FLAG} \
	--metadata-dir "${metadata_dir}" \
	--files-dir "${unencrypted_files_dir}" \
	--output-encrypted-files-dir "${output_encrypted_files_dir}" \
	--logs-dir "${progress_logs_dir}" \
	--force \
	>>"$log_stdout" 2>>"$log_stderr"
