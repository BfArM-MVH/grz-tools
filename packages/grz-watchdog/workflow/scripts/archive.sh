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
submission_dir="${snakemake_input[data]}"

CONSENT=$(cat "${snakemake_input[consent_flag]}")
if [[ "$CONSENT" == "true" ]]; then
	CONFIG_FILE="${snakemake_input[consented_config_path]}"
else
	CONFIG_FILE="${snakemake_input[nonconsented_config_path]}"
fi

echo "Consent: $CONSENT. Using config file for archiving: $CONFIG_FILE" >>"$log_stdout" 2>>"$log_stderr"

echo "Redacting sensitive data from logs before archiving..." >>"$log_stdout" 2>>"$log_stderr"
metadata_file="${submission_dir}/metadata/metadata.json"
logs_dir="${submission_dir}/logs"

if [ -f "$metadata_file" ] && [ -d "$logs_dir" ]; then
	tanG=$(jq --raw-output '.submission.tanG' "$metadata_file" || true)
	localCaseId=$(jq --raw-output '.submission.localCaseId' "$metadata_file" || true)

	if [ -n "$tanG" ] && [ "$tanG" != "null" ]; then
		echo "Redacting tanG..." >>"$log_stdout" 2>>"$log_stderr"
		rg --files-with-matches --fixed-strings "$tanG" "$logs_dir" | xargs --no-run-if-empty sed -i "s/$tanG/0000000000000000000000000000000000000000000000000000000000000000/g" || true
	fi

	if [ -n "$localCaseId" ] && [ "$localCaseId" != "null" ]; then
		echo "Redacting localCaseId..." >>"$log_stdout" 2>>"$log_stderr"
		rg --files-with-matches --fixed-strings "$localCaseId" "$logs_dir" | xargs --no-run-if-empty sed -i "s/$localCaseId/REDACTED_LOCAL_CASE_ID/g" || true
	fi
	echo "Redaction complete." >>"$log_stdout" 2>>"$log_stderr"
else
	echo "[WARNING] metadata.json or logs directory not found. Skipping log redaction." >>"$log_stdout" 2>>"$log_stderr"
fi

grzctl db --config-file "${db_config}" submission update --ignore-error-state "${submission_id}" archiving >>"$log_stdout" 2>>"$log_stderr"

grzctl archive \
	--config-file "$CONFIG_FILE" \
	--submission-dir "${submission_dir}" \
	>>"$log_stdout" 2>>"$log_stderr"

grzctl db --config-file "${db_config}" submission update --ignore-error-state "${submission_id}" archived >>"$log_stdout" 2>>"$log_stderr"
