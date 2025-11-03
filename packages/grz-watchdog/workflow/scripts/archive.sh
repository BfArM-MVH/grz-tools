#!/usr/bin/env bash
set -euo pipefail

submission_id="${snakemake_wildcards[submission_id]}"
db_config="${snakemake_input[db_config_path]}"
log_stdout="${snakemake_log[stdout]}"
log_stderr="${snakemake_log[stderr]}"

metadata_file_path="${snakemake_input[metadata]}"
re_encrypted_data_dir="${snakemake_input[re_encrypted_data]}"
consent_flag_path="${snakemake_input[consent_flag]}"
consented_config_path="${snakemake_input[consented_config_path]}"
nonconsented_config_path="${snakemake_input[nonconsented_config_path]}"

read -r -a logs_to_archive <<<"${snakemake_input[logs_to_archive]}"

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

grzctl db --config-file "${db_config}" submission update --ignore-error-state "${submission_id}" archiving >"$log_stdout" 2>"$log_stderr"

staging_dir=$(mktemp -d -p "$(dirname ${snakemake_output[marker]})" "archive_staging_XXXXXX")
trap 'rm -rf "$staging_dir"' EXIT
echo "Created temporary archive staging directory: $staging_dir" >>"$log_stdout"

echo "Assembling archive package..." >>"$log_stdout"
ln -s "$(realpath "${re_encrypted_data_dir}")/encrypted_files" "${staging_dir}/encrypted_files"
mkdir "${staging_dir}/metadata"
cp "${metadata_file_path}" "${staging_dir}/metadata/metadata.json"
mkdir "${staging_dir}/logs"
cp "${logs_to_archive[@]}" "${staging_dir}/logs/"

echo "Redacting sensitive data from staged logs..." >>"$log_stdout"
staged_metadata_file="${staging_dir}/metadata/metadata.json"
staged_logs_dir="${staging_dir}/logs"

tanG=$(jq --raw-output '.submission.tanG' "$staged_metadata_file" || true)
localCaseId=$(jq --raw-output '.submission.localCaseId' "$staged_metadata_file" || true)

if [ -n "$tanG" ] && [ "$tanG" != "null" ]; then
	rg --files-with-matches --fixed-strings "$tanG" "$staged_logs_dir" | xargs --no-run-if-empty sed -i "s/$tanG/$(printf '0'%.0s {1..64})/g" || true
fi
if [ -n "$localCaseId" ] && [ "$localCaseId" != "null" ]; then
	rg --files-with-matches --fixed-strings "$localCaseId" "$staged_logs_dir" | xargs --no-run-if-empty sed -i "s/$localCaseId/REDACTED_LOCAL_CASE_ID/g" || true
fi
echo "Redaction complete." >>"$log_stdout"

CONSENT=$(cat "${consent_flag_path}")
if [[ "$CONSENT" == "true" ]]; then
	CONFIG_FILE="${consented_config_path}"
else
	CONFIG_FILE="${nonconsented_config_path}"
fi

echo "Consent: $CONSENT. Using config file for archiving: $CONFIG_FILE" >>"$log_stdout"

grzctl archive \
	--config-file "$CONFIG_FILE" \
	--submission-dir "${staging_dir}" \
	>>"$log_stdout" 2>>"$log_stderr"

grzctl db --config-file "${db_config}" submission update --ignore-error-state "${submission_id}" archived >>"$log_stdout" 2>>"${log_stderr}"
