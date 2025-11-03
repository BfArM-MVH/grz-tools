#!/usr/bin/env bash
set -euo pipefail

_error_handler() {
	local exit_code="$1"
	local line_no="$2"
	local command="$3"

	local error_message="[ERROR] Script '$0' failed on line $line_no with exit code $exit_code while executing command: $command"
	echo "$error_message" >&2
	echo "$error_message" >>"${log_stderr}"

	popd # needed so relative paths specified in db_config can be resolved correctly, since nextflow is called from within $launch_dir
	grzctl db --config-file "${db_config}" submission update --ignore-error-state "${submission_id}" error >>"${log_stdout}" 2>>"${log_stderr}"
}

trap '_error_handler $? $LINENO "$BASH_COMMAND"' ERR

submission_id="${snakemake_wildcards[submission_id]}"
db_config="${snakemake_input[db_config_path]}"
log_stdout=$(realpath "${snakemake_log[stdout]}")
log_stderr=$(realpath "${snakemake_log[stderr]}")

launch_dir=$(realpath "${snakemake_input[launch_dir]}")
work_dir=$(realpath "${snakemake_output[work_dir]}")
out_dir=$(realpath "${snakemake_output[out_dir]}")

pipeline=$(realpath "${snakemake_input[pipeline]}")
reference_path=$(realpath "${snakemake_input[reference_path]}")

decrypted_data_dir_abs=$(realpath "${snakemake_input[submission_basepath]}")
metadata_file_abs=$(realpath "${snakemake_input[metadata]}")

configs="${snakemake_params[configs]}"
profiles="${snakemake_params[profiles]}"

extra="${snakemake_params[extra]}"

grzctl db --config-file "${db_config}" submission update --ignore-error-state "${submission_id}" qcing >>"$log_stdout" 2>>"$log_stderr"

mkdir -p "${work_dir}"
mkdir -p "${out_dir}"
pushd "${launch_dir}"

qc_submission_stage_dir="qc_stage_${submission_id}"
mkdir -p "$qc_submission_stage_dir"
trap 'rm -rf "$qc_submission_stage_dir"' EXIT SIGHUP SIGINT SIGTERM

ln -sfn "${decrypted_data_dir_abs}/files" "${qc_submission_stage_dir}/files"
mkdir -p "${qc_submission_stage_dir}/metadata"
ln -sfn "${metadata_file_abs}" "${qc_submission_stage_dir}/metadata/metadata.json"

submission_basepath_for_nextflow="$qc_submission_stage_dir"

nextflow run "${pipeline}" \
	${configs} \
	-profile ${profiles} \
	--outdir "${out_dir}" \
	--reference_path "${reference_path}" \
	--submission_basepath "${submission_basepath_for_nextflow}" \
	-work-dir "${work_dir}" \
	${extra} \
	>>"$log_stdout" 2>>"$log_stderr"

popd
grzctl db --config-file "${db_config}" submission update --ignore-error-state "${submission_id}" qced >>"$log_stdout" 2>>"$log_stderr"
