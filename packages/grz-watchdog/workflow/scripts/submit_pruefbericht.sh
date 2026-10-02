#!/usr/bin/env bash
set -euo pipefail

submission_id="${snakemake_wildcards[submission_id]}"
grzctl_config="${snakemake_input[grzctl_config_path]}"
log_stdout="${snakemake_log[stdout]}"
log_stderr="${snakemake_log[stderr]}"

# if snakemake_params[custom_ca_cert] exists, prepend REQUESTS_CA_BUNDLE=…
if [ -n "${snakemake_params[custom_ca_cert]:-}" ]; then
	export REQUESTS_CA_BUNDLE="${snakemake_params[custom_ca_cert]}"
fi

# grzctl pruefbericht submit handles DB state transitions (REPORTING → REPORTED) via DbContext.
grzctl --config "${grzctl_config}" pruefbericht \
	submit \
	--submission-id "${submission_id}" \
	--pruefbericht-file "${snakemake_input[pruefbericht]}" \
	--print-token \
	>"${snakemake_output[answer]}" 2>>"$log_stderr"
