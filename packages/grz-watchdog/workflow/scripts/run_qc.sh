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

launch_dir="${snakemake_input[launch_dir]}"

work_dir="${snakemake_output[work_dir]}"

out_dir="${snakemake_params[out_dir]}"
absolute_pipeline_path="${snakemake_params[absolute_pipeline_path]}"
absolute_submission_basepath="${snakemake_params[absolute_submission_basepath]}"
reference_path="${snakemake_params[reference_path]}"
configs="${snakemake_params[configs]}"
profiles="${snakemake_params[profiles]}"

extra="${snakemake_params[extra]}"

grzctl db --config-file "${db_config}" submission update --ignore-error-state "${submission_id}" qcing >>"$log_stdout" 2>>"$log_stderr"

mkdir -p "${work_dir}"
mkdir -p "${out_dir}"
cd "${launch_dir}"
nextflow run "${absolute_pipeline_path}" \
	"${configs}" \
	-profile "${profiles}" \
	--outdir "${out_dir}" \
	--reference_path "${reference_path}" \
	--submission_basepath "${absolute_submission_basepath}" \
	"${extra}" \
	>>"$log_stdout" 2>>"$log_stderr"

grzctl db --config-file "${db_config}" submission update --ignore-error-state "${submission_id}" qced >>"$log_stdout" 2>>"$log_stderr"