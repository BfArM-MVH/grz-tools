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

if [ -n "${snakemake_params[s3_access_key]:-}" ]; then
	export GRZ_S3__ACCESS_KEY="${snakemake_params[s3_access_key]}"
fi

if [ -n "${snakemake_params[s3_secret]:-}" ]; then
	export GRZ_S3__SECRET="${snakemake_params[s3_secret]}"
fi

if [ -n "${snakemake_params[db_author_key_passphrase]:-}" ]; then
	export GRZ_DB__AUTHOR__PRIVATE_KEY_PASSPHRASE="${snakemake_params[db_author_key_passphrase]}"
fi

submission_id="${snakemake_wildcards[submission_id]}"
db_config="${snakemake_input[db_config_path]}"
log_stdout="${snakemake_log[stdout]}"
log_stderr="${snakemake_log[stderr]}"

grzctl db --config-file "${db_config}" submission update --ignore-error-state "${submission_id}" downloading >"$log_stdout" 2>"$log_stderr"

grzctl download \
	--submission-id "${submission_id}" \
	--output-dir "${snakemake_output[data]}" \
	--config-file "${snakemake_input[inbox_config_path]}" \
	>>"$log_stdout" 2>>"$log_stderr"

grzctl db --config-file "${db_config}" submission update --ignore-error-state "${submission_id}" downloaded >>"$log_stdout" 2>>"$log_stderr"

grzctl db --config-file "${db_config}" submission populate --no-confirm "${submission_id}" "${snakemake_output[data]}/metadata/metadata.json" >>"$log_stdout" 2>>"$log_stderr"
