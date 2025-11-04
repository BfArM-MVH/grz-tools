#!/usr/bin/env bash
set -euo pipefail

_error_handler() {
	local exit_code="$1"
	local line_no="$2"
	local command="$3"

	local error_message="[ERROR] Script '$0' failed on line $line_no with exit code $exit_code while executing command: $command"
	echo "$error_message" >&2
	echo "$error_message" >>"${log_stderr}"

	grzctl db --config-file "${db_config}" submission update --ignore-error-state "${submission_id}" error >>"${log_stdout}" 2>>"${log_stderr}"
}

trap '_error_handler $? $LINENO "$BASH_COMMAND"' ERR

if [ -n "${snakemake_params[db_author_key_passphrase]:-}" ]; then
	export GRZ_DB__AUTHOR__PRIVATE_KEY_PASSPHRASE="${snakemake_params[db_author_key_passphrase]}"
fi

if [ -n "${snakemake_params[grz_private_key_passphrase]:-}" ]; then
	export C4GH_PASSPHRASE="${snakemake_params[grz_private_key_passphrase]}"
fi

submission_id="${snakemake_wildcards[submission_id]}"
db_config="${snakemake_input[db_config_path]}"
log_stdout="${snakemake_log[stdout]}"
log_stderr="${snakemake_log[stderr]}"
inbox_config_path="${snakemake_input[inbox_config_path]}"

input_metadata_file="${snakemake_input[metadata]}"
metadata_dir="$(dirname "$metadata_file_path")"
input_encrypted_files_dir="${snakemake_input[encrypted_files_dir]}"
output_files_dir="${snakemake_output[files_dir]}"
output_logs_dir="${snakemake_output[logs_dir]}"
progress_logs_dir="$(dirname "${snakemake_output[progress_log]}")"

grzctl db --config-file "${db_config}" submission update --ignore-error-state "${submission_id}" decrypting >"$log_stdout" 2>"$log_stderr"

grzctl decrypt \
	--config-file "${inbox_config_path}" \
	--metadata-dir "${metadata_dir}" \
	--encrypted-files-dir "${input_encrypted_files_dir}" \
	--files-dir "${output_files_dir}" \
	--logs-dir "${progress_logs_dir}" \
	>>"$log_stdout" 2>>"$log_stderr"

grzctl db --config-file "${db_config}" submission update --ignore-error-state "${submission_id}" decrypted >>"${log_stdout}" 2>>"${log_stderr}"
