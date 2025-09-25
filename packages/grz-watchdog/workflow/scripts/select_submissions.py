import json
import subprocess
from datetime import datetime
from contextlib import redirect_stdout, redirect_stderr


def select_submissions(submissions_path, db_config_path, limit):
    result = subprocess.run(
        ["grzctl", "db", "--config-file", db_config_path, "list", "--json"],
        check=True,
        capture_output=True,
        text=True,
    )
    db_output = result.stdout
    if not db_output:
        print("No submissions found in database. Nothing to do.")
        return []
    db_states = {}
    for db_submission_entry in json.loads(db_output):
        if latest_state := db_submission_entry.get("latest_state"):
            db_states[db_submission_entry["id"]] = {
                "state": latest_state["state"].casefold(),
                "timestamp": latest_state["timestamp"],
            }

    with open(submissions_path, "r") as f:
        synced_submissions = json.load(f)

    pending_submissions = []
    for submission in synced_submissions:
        submission_id = submission["submission_id"]
        latest_db_state = db_states.get(submission_id)

        if latest_db_state and latest_db_state["state"] == "uploaded":
            submission["timestamp"] = latest_db_state["timestamp"]
            pending_submissions.append(submission)

    pending_submissions.sort(key=lambda s: datetime.fromisoformat(s["timestamp"]))
    if limit is not None:
        pending_submissions = pending_submissions[:limit]

    return pending_submissions


submissions_path = snakemake.input.submissions
db_config_path = snakemake.input.db_config_path
output_batch_file = snakemake.output.submissions_batch
limit = snakemake.params.batch_limit
stdout_log_path = snakemake.log.stdout
stderr_log_path = snakemake.log.stderr

with (
    redirect_stdout(open(stdout_log_path, "w")),
    redirect_stderr(open(stderr_log_path, "w")),
):
    selected_submissions = select_submissions(submissions_path, db_config_path, limit)
    print(f"Selected {len(selected_submissions)} submissions for processing.")
    with open(output_batch_file, "w") as f:
        for submission in selected_submissions:
            origin = submission["origin"]
            # TODO implement selection strategy for with_qc/without_qc
            target_path = (
                f"results/{origin['submitter_id']}/{origin['inbox']}/{submission['submission_id']}/processed/without_qc"
            )
            f.write(f"{target_path}\n")
