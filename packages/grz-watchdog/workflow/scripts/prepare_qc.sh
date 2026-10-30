#!/usr/bin/env bash
set -euo pipefail

_error_handler() {
	local exit_code="$1"
	local line_no="$2"
	local command="$3"

	local error_message="[ERROR] Script '$0' failed on line $line_no with exit code $exit_code while executing command: $command"
	echo "$error_message" >&2
	echo "$error_message" >>"${log_stderr}"
}

trap '_error_handler $? $LINENO "$BASH_COMMAND"' ERR

log_stdout=$(realpath "${snakemake_log[stdout]}")
log_stderr=$(realpath "${snakemake_log[stderr]}")

launch_dir=$(realpath "${snakemake_output[launch_dir]}")
work_dir=$(realpath "${snakemake_output[work_dir]}")

pipeline=$(realpath "${snakemake_input[pipeline]}")
reference_path=$(realpath "${snakemake_output[references_dir]}")
out_dir="${reference_path}"

configs="${snakemake_params[configs]}" # these already are absolute paths
profiles="${snakemake_params[profiles]}"

extra="${snakemake_params[extra]}"

mkdir -p "${work_dir}"
mkdir -p "${out_dir}"
mkdir -p "${launch_dir}"
pushd "${launch_dir}"

nextflow run "${pipeline}" \
	${configs} \
	-profile ${profiles} \
	--outdir "${out_dir}" \
	--reference_path "${reference_path}" \
	-work-dir "${work_dir}" \
	${extra} \
	>>"$log_stdout" 2>>"$log_stderr"

popd
