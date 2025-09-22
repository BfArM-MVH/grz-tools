#!/usr/bin/env bash
set -euo pipefail

log_stdout="${snakemake_log[stdout]}"
log_stderr="${snakemake_log[stderr]}"

(
	IS_VALID=$(cat "${snakemake_input[validation_flag]}")
	if [[ "$IS_VALID" == "true" ]]; then
		pruefbericht_params=""
	else
		pruefbericht_params="--fail"
	fi

	grzctl pruefbericht \
		--config-file "${snakemake_input[pruefbericht_config_path]}" \
		--submission-dir "${snakemake_input[data]}" \
		${pruefbericht_params} \
		--dry-run --json >"${snakemake_output[pruefbericht]}"
) >"$log_stdout" 2>"$log_stderr"
