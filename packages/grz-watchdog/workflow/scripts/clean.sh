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

artifacts_to_remove=(
	"${snakemake_input[downloaded_data]}"
	"${snakemake_input[decrypted_data]}"
	"${snakemake_input[re_encrypted_data]}"
	"${snakemake_input[metadata_file]}"
	"${snakemake_input[validation_flag]}"
	"${snakemake_input[validation_errors]}"
	"${snakemake_input[consent_flag]}"
)

read -r -a progress_log_files <<<"${snakemake_input[progress_logs]}"

if [ ${#progress_log_files[@]} -gt 0 ]; then
	progress_logs_dir="$(dirname "${progress_log_files[0]}")"
	artifacts_to_remove+=("$progress_logs_dir")
fi

# If mode is 'none', do nothing and exit successfully.
if [[ "$mode" == "none" ]]; then
	echo "Auto-cleanup mode is 'none'. No action taken." >>"$log_stdout"
	echo 'true' >"${snakemake_output[clean_results]}"
	exit 0
fi

grzctl db --config-file "${db_config}" submission update --ignore-error-state "${submission_id}" cleaning >"$log_stdout" 2>"$log_stderr"

echo "Auto-cleanup mode: '${mode}'" >>"$log_stdout"

if [[ "$mode" == "none" ]]; then
	echo "Auto-cleanup mode is 'none'. No local storage action taken." >>"$log_stdout"
else
	if [[ "$mode" == "storage" || "$mode" == "inbox+storage" ]]; then
		echo "Cleaning all local submission artifacts..." >>"$log_stdout"
		for item in "${artifacts_to_remove[@]}"; do
			if [[ -e "$item" ]]; then
				echo "Removing: $item" >>"$log_stdout"
				rm -rf "$item"
			else
				echo "[WARNING] Not found, skipping: $item" >>"$log_stdout"
			fi
		done
		echo "Local storage successfully removed." >>"$log_stdout"
	fi
fi

if [[ "$mode" == "inbox" || "$mode" == "inbox+storage" ]]; then
	echo "Cleaning S3 inbox..." >>"$log_stdout"
	grzctl db --config-file "${db_config}" submission update --ignore-error-state "${submission_id}" cleaning >"$log_stdout" 2>"$log_stderr"
	grzctl clean --config-file "${inbox_config}" --submission-id "${submission_id}" --yes-i-really-mean-it >>"$log_stdout" 2>>"${log_stderr}"
fi

details=""
if [[ "$mode" == "inbox" || "$mode" == "inbox+storage" ]]; then
	details='"inbox"'
fi
if [[ "$mode" == "storage" || "$mode" == "inbox+storage" ]]; then
	if [[ -n "$details" ]]; then details+=", "; fi
	details+='"local_storage"'
fi
json_data="{\"targets\": [${details}]}"

if [[ "$mode" != "none" ]]; then
	grzctl db --config-file "${db_config}" submission update --ignore-error-state "${submission_id}" cleaned --data "$json_data" >>"$log_stdout" 2>>"$log_stderr"
fi

echo 'true' >"${snakemake_output[clean_results]}"
