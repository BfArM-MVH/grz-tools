#!/usr/bin/env bash
set -euo pipefail

submission_id="${snakemake_wildcards[submission_id]}"
grzctl_config="${snakemake_input[grzctl_config_path]}"

(
	echo "Submission ${submission_id} failed validation."
	# grzctl validate has recorded the failure with its reason, but a cleanup records later states.
	# Record the failure again, so that the ERROR state keeps the submission out of later batches.
	states=$(grzctl --config "${grzctl_config}" db submission show --json "${submission_id}" | jq --compact-output '.states')
	if [[ $(jq --raw-output 'last.state' <<<"$states") != "Error" ]]; then
		reason=$(jq --raw-output 'map(select(.state == "Error")) | last | .failure_reason // "unknown"' <<<"$states")
		grzctl --config "${grzctl_config}" db submission update --ignore-error-state "${submission_id}" error --failure-reason "$reason"
	fi
	echo "Submission ${submission_id} processing finished due to validation failure." >"${snakemake_output[target]}"
) >"${snakemake_log[stdout]}" 2>"${snakemake_log[stderr]}"
