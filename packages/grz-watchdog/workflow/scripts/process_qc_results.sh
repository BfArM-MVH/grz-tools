#!/usr/bin/env bash
set -euo pipefail

_error_handler() {
	local exit_code="$1"
	local line_no="$2"
	local command="$3"

	local error_message="[ERROR] Script '$0' failed on line $line_no with exit code $exit_code while executing command: $command"
	echo "$error_message" >&2
	echo "$error_message" >>"${log_stderr}"

	grzctl --config "${grzctl_config}" db submission update --ignore-error-state "${submission_id}" error --failure-reason detailed_qc_error >>"${log_stdout}" 2>>"${log_stderr}"
}

trap '_error_handler $? $LINENO "$BASH_COMMAND"' ERR

submission_id="${snakemake_wildcards[submission_id]}"
grzctl_config="${snakemake_input[grzctl_config_path]}"
report_csv="${snakemake_params[report_csv]}"
qc_workflow_version="${snakemake_params[qc_workflow_version]}"
log_stdout="${snakemake_log[stdout]}"
log_stderr="${snakemake_log[stderr]}"

# Since GRZ_QC_Workflow 4.0.0, a deviation from the values that the LE provided is reported as DEVIATION.
# A deviation does not fail the QC, so the QC passes if every lab datum of the index donor has PASS or DEVIATION.
qc_status=$(awk -F, '$2 == "index"' <"${report_csv}" | python3 -c 'import csv,sys; s={row[6] for row in csv.reader(sys.stdin)}; print(str(bool(s) and s <= {"PASS", "DEVIATION"}).lower())')

grzctl --config "${grzctl_config}" db submission modify "${submission_id}" detailed_qc_passed "${qc_status}" >"$log_stdout" 2>"$log_stderr"
grzctl --config "${grzctl_config}" db submission populate-qc --no-confirm "${submission_id}" "${report_csv}" >>"$log_stdout" 2>>"$log_stderr"
