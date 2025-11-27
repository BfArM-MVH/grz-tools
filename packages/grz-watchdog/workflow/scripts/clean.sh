#!/usr/bin/env bash
set -euo pipefail

_error_handler() {
	local exit_code="$1"
	local line_no="$2"
	local command="$3"

	local error_message="[ERROR] Script '$0' failed on line $line_no with exit code $exit_code while executing command: $command"
	echo "$error_message" >&2
	echo "$error_message" >>"${log_stderr}"

	# If the script fails, we should still record an error state in the DB
	grzctl db --config-file "${db_config}" submission update "${submission_id}" error --data '{"reason": "cleanup script failed"}' >>"${log_stdout}" 2>>"${log_stderr}"
}

trap '_error_handler $? $LINENO "$BASH_COMMAND"' ERR

submission_id="${snakemake_wildcards[submission_id]}"
inbox_config="${snakemake_input[inbox_config_path]}"
db_config="${snakemake_input[db_config_path]}"
log_stdout="${snakemake_log[stdout]}"
log_stderr="${snakemake_log[stderr]}"
mode="${snakemake_params[mode]}"

# If mode is 'none', do nothing and exit successfully.
if [[ "$mode" == "none" ]]; then
	echo "Auto-cleanup mode is 'none'. No action taken." >>"$log_stdout"
	echo 'true' >"${snakemake_output[clean_results]}"
	exit 0
fi

echo "Auto-cleanup mode: '${mode}'" >>"$log_stdout"

details=""
if [[ "$mode" == "inbox" ]]; then
	echo "Cleaning S3 inbox..." >>"$log_stdout"
	grzctl db --config-file "${db_config}" submission update --ignore-error-state "${submission_id}" cleaning >"$log_stdout" 2>"$log_stderr"
	grzctl clean --config-file "${inbox_config}" --submission-id "${submission_id}" --yes-i-really-mean-it >>"$log_stdout" 2>>"${log_stderr}"
	details='"inbox"'
fi

json_data="{\"targets\": [${details}]}"

if [[ "$mode" != "none" ]]; then
	grzctl db --config-file "${db_config}" submission update --ignore-error-state "${submission_id}" cleaned --data "$json_data" >>"$log_stdout" 2>>"$log_stderr"
fi

echo 'true' >"${snakemake_output[clean_results]}"
