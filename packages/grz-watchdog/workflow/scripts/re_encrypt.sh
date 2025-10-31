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

submission_id="${snakemake_wildcards[submission_id]}"
db_config="${snakemake_input[db_config_path]}"
log_stdout="${snakemake_log[stdout]}"
log_stderr="${snakemake_log[stderr]}"

CONSENT=$(cat "${snakemake_input[consent_flag]}")
if [[ "$CONSENT" == "true" ]]; then
	CONFIG_FILE="${snakemake_input[consented_config_path]}"
else
	CONFIG_FILE="${snakemake_input[nonconsented_config_path]}"
fi

echo "Consent: $CONSENT. Using config file: $CONFIG_FILE" >"$log_stdout" 2>"$log_stderr"

grzctl db --config-file "${db_config}" submission update --ignore-error-state "${submission_id}" encrypting >>"$log_stdout" 2>>"$log_stderr"

grzctl encrypt \
	--config-file "$CONFIG_FILE" \
	--submission-dir "${snakemake_input[data]}" \
	--force \
	>>"$log_stdout" 2>>"$log_stderr"

grzctl db --config-file "${db_config}" submission update --ignore-error-state "${submission_id}" encrypted >>"$log_stdout" 2>>"$log_stderr"
