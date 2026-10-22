#!/usr/bin/env bash
set -euo pipefail

submission_id="${snakemake_wildcards[submission_id]}"
inbox="${snakemake_wildcards[inbox]}"
grzctl_config="${snakemake_input[grzctl_config_path]}"
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

if [[ "$mode" == "inbox" ]]; then
	echo "Cleaning S3 inbox..." >>"$log_stdout"
	# grzctl clean handles DB state transitions (CLEANING → CLEANED, or ERROR on failure) via DbContext.
	grzctl --config "${grzctl_config}" clean --inbox "${inbox}" --submission-id "${submission_id}" --yes-i-really-mean-it >>"$log_stdout" 2>>"${log_stderr}"
fi

echo 'true' >"${snakemake_output[clean_results]}"
