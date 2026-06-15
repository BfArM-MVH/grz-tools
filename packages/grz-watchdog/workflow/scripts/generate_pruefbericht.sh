#!/usr/bin/env bash
set -euo pipefail

log_stdout="${snakemake_log[stdout]}"
log_stderr="${snakemake_log[stderr]}"
db_config="${snakemake_input[db_config_path]}"
submission_id="${snakemake_wildcards[submission_id]}"

(
	IS_VALID=$(cat "${snakemake_input[validation_flag]}")
	if [[ "$IS_VALID" == "true" ]]; then
		pruefbericht_params=""
	else
		pruefbericht_params="--fail"
	fi

	grzctl --config-file "${db_config}" pruefbericht \
		generate \
		from-database \
		--submission-id "${submission_id}" \
		${pruefbericht_params} \
		>"${snakemake_output[pruefbericht]}"
) >"$log_stdout" 2>"$log_stderr"
