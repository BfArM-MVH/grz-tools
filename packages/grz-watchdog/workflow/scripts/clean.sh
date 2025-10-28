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
local_data_dir="${snakemake_input[data]}"
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

grzctl db --config-file "${db_config}" submission update --ignore-error-state "${submission_id}" cleaning >"$log_stdout" 2>"$log_stderr"

echo "Auto-cleanup mode: '${mode}'" >>"$log_stdout"

case "$mode" in
"inbox")
	echo "Cleaning S3 inbox..." >>"$log_stdout"
	grzctl clean --config-file "${inbox_config}" --submission-id "${submission_id}" --yes-i-really-mean-it >>"$log_stdout" 2>>"$log_stderr"
	;;

"storage")
	echo "Cleaning local storage at: ${local_data_dir}" >>"$log_stdout"
	if [[ -d "$local_data_dir" ]]; then
		rm -rf "$local_data_dir"
		echo "Local storage successfully removed." >>"$log_stdout"
	else
		echo "[WARNING] Local storage directory not found, skipping." >>"$log_stdout"
	fi
	;;

"inbox+storage")
	echo "Cleaning S3 inbox..." >>"$log_stdout"
	grzctl clean --config-file "${inbox_config}" --submission-id "${submission_id}" --yes-i-really-mean-it >>"$log_stdout" 2>>"$log_stderr"

	echo "Cleaning local storage at: ${local_data_dir}" >>"$log_stdout"
	if [[ -d "$local_data_dir" ]]; then
		rm -rf "$local_data_dir"
		echo "Local storage successfully removed." >>"$log_stdout"
	else
		echo "[WARNING] Local storage directory not found, skipping." >>"$log_stdout"
	fi
	;;

*)
	echo "[ERROR] Unknown auto-cleanup mode: '${mode}'" >>"$log_stderr"
	exit 1
	;;
esac

details=""
if [[ "$mode" == "inbox" || "$mode" == "inbox+storage" ]]; then
	details='"inbox"'
fi
if [[ "$mode" == "storage" || "$mode" == "inbox+storage" ]]; then
	if [[ -n "$details" ]]; then details+=", "; fi
	details+='"local_storage"'
fi
json_data="{\"targets\": [${details}]}"

grzctl db --config-file "${db_config}" submission update --ignore-error-state "${submission_id}" cleaned --data "$json_data" >>"$log_stdout" 2>>"$log_stderr"

echo 'true' >"${snakemake_output[clean_results]}"
