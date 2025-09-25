#!/usr/bin/env bash
set -euo pipefail

_error_handler() {
	local exit_code="$1"
	local line_no="$2"
	local command="$3"

	local error_message="[ERROR] Script '$0' failed on line $line_no with exit code $exit_code while executing command: $command"
	echo "$error_message" >&2
	echo "$error_message" >>"${log_stderr}"

	grzctl db --config-file "${db_config}" submission update "${submission_id}" error >>"${log_stdout}" 2>>"${log_stderr}"
}

trap '_error_handler $? $LINENO "$BASH_COMMAND"' ERR

submission_id="${snakemake_wildcards[submission_id]}"
inbox_config="${snakemake_input[inbox_config_path]}"
db_config="${snakemake_input[db_config_path]}"
log_stdout="${snakemake_log[stdout]}"
log_stderr="${snakemake_log[stderr]}"
mode="${snakemake_params[mode]}"  # TODO actually respect the mode

grzctl db --config-file "$db_config" submission update "$submission_id" cleaning >"$log_stdout" 2>"$log_stderr"
grzctl clean --config-file "$inbox_config" --submission-id "$submission_id" --yes-i-really-mean-it >>"$log_stdout" 2>>"$log_stderr"
echo 'true' >"${snakemake_output[clean_results]}"
grzctl db --config-file "$db_config" submission update "$submission_id" cleaned >>"$log_stdout" 2>>"$log_stderr"
