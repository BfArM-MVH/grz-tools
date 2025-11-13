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
		generate \
		from-metadata "${snakemake_input[metadata]}" \
		${pruefbericht_params} \
		>"${snakemake_output[pruefbericht]}"
) >"$log_stdout" 2>"$log_stderr"
