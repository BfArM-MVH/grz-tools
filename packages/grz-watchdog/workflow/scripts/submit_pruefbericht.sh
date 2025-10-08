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

# if snakemake_params[custom_ca_cert] exists, prepend REQUESTS_CA_BUNDLE=…
if [ -n "${snakemake_params[custom_ca_cert]:-}" ]; then
	export REQUESTS_CA_BUNDLE="${snakemake_params[custom_ca_cert]}"
fi

grzctl pruefbericht \
	submit \
	--config-file "${snakemake_input[pruefbericht_config_path]}" \
	--pruefbericht-file "${snakemake_input[pruefbericht]}" \
	--print-token \
	>"${snakemake_output[answer]}" 2>>"$log_stderr"

grzctl db --config-file "${db_config}" submission update --ignore-error-state "${submission_id}" reported >>"$log_stdout" 2>>"$log_stderr"
