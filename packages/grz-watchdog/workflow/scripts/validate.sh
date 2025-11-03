#!/usr/bin/env bash
set -euo pipefail

submission_id="${snakemake_wildcards[submission_id]}"
db_config="${snakemake_input[db_config_path]}"
inbox_config="${snakemake_input[inbox_config_path]}"

metadata_file_path="${snakemake_input[metadata]}"
metadata_dir="$(dirname "$metadata_file_path")"
decrypted_data_dir="${snakemake_input[data]}"
progress_logs_dir="$(dirname "${snakemake_output[checksum_log]}")"
mkdir -p "${progress_logs_dir}"

validation_flag="${snakemake_output[validation_flag]}"
validation_errors="${snakemake_output[validation_errors]}"
log_stdout="${snakemake_log[stdout]}"
log_stderr="${snakemake_log[stderr]}"

grzctl db --config-file "${db_config}" submission update --ignore-error-state "${submission_id}" validating >"$log_stdout" 2>"$log_stderr"

# We expect `grzctl validate` to return a non-zero code on validation failure,
# which is not a script error. So we handle its exit code manually instead of relying on `set -e`.
if grzctl validate \
	--config-file "${inbox_config}" \
	--metadata-dir "${metadata_dir}" \
	--files-dir "${decrypted_data_dir}/files" \
	--logs-dir "${progress_logs_dir}" \
	>>"$log_stdout" 2>"$validation_errors"; then
	echo "true" >"$validation_flag"
	grzctl db --config-file "${db_config}" submission modify "${submission_id}" basic_qc_passed yes
else
	# Failure: Validation found errors. This is an expected outcome.
	# The errors are already captured in the validation_errors file.
	echo "false" >"$validation_flag"
	grzctl db --config-file "${db_config}" submission modify "${submission_id}" basic_qc_passed no
fi

grzctl db --config-file "${db_config}" submission update --ignore-error-state "${submission_id}" validated >>"$log_stdout" 2>>"$log_stderr"
