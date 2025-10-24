#!/usr/bin/env bash
set -euo pipefail

submission_id="${snakemake_wildcards[submission_id]}"
db_config="${snakemake_input[db_config_path]}"
inbox_config="${snakemake_input[inbox_config_path]}"
submission_dir="${snakemake_input[data]}"
validation_flag="${snakemake_output[validation_flag]}"
validation_errors="${snakemake_output[validation_errors]}"
log_stdout="${snakemake_log[stdout]}"
log_stderr="${snakemake_log[stderr]}"

grzctl db --config-file "${db_config}" submission update --ignore-error-state "${submission_id}" validating >"$log_stdout" 2>"$log_stderr"

# We expect `grzctl validate` to return a non-zero code on validation failure,
# which is not a script error. So we handle its exit code manually instead of relying on `set -e`.
if grzctl validate \
	--config-file "${inbox_config}" \
	--submission-dir "$submission_dir" \
	>>"$log_stdout" 2>"$validation_errors"; then
	echo "true" >"$validation_flag"
else
	# Failure: Validation found errors. This is an expected outcome.
	# The errors are already captured in the validation_errors file.
	echo "false" >"$validation_flag"
fi

grzctl db --config-file "${db_config}" submission update --ignore-error-state "${submission_id}" validated >>"$log_stdout" 2>>"$log_stderr"
