import json
import os
import sys
from contextlib import redirect_stderr, redirect_stdout
from datetime import datetime

sys.path.append(os.path.dirname(__file__))
import shared


def select_submissions(inbox_scans, db_config_path, limit):
    all_s3_submissions = []
    for inbox_scan in inbox_scans:
        with open(inbox_scan) as f:
            all_s3_submissions.extend(json.load(f))

    db_states = shared.get_db_states(db_config_path)

    pending_submissions = []
    for submission in all_s3_submissions:
        submission_id = submission["submission_id"]
        s3_state = submission.get("state")
        latest_db_state = db_states.get(submission_id)

        if latest_db_state and latest_db_state.get("timestamp"):
            submission["timestamp"] = latest_db_state["timestamp"]
        else:
            s3_timestamp_str = submission.get("latest_modification", "1970-01-01T00:00:00Z")
            dt_object = datetime.fromisoformat(s3_timestamp_str.replace("Z", "+00:00"))
            submission["timestamp"] = dt_object.isoformat()

        if (
            latest_db_state
            and (latest_db_state["state"] == "uploaded" or latest_db_state["state"] in shared.CONTINUABLE_STATES)
        ) or (not latest_db_state and s3_state == "complete"):
            pending_submissions.append(submission)

    pending_submissions.sort(key=lambda s: datetime.fromisoformat(s["timestamp"]))
    if limit is not None:
        pending_submissions = pending_submissions[:limit]

    return pending_submissions


inbox_scans = snakemake.input.scans
db_config_path = snakemake.input.db_config_path
output_batch_file = snakemake.output.submissions_batch
limit = snakemake.params.batch_limit
stdout_log_path = snakemake.log.stdout
stderr_log_path = snakemake.log.stderr

with redirect_stdout(open(stdout_log_path, "w")), redirect_stderr(open(stderr_log_path, "w")):
    selected_submissions = select_submissions(inbox_scans, db_config_path, limit)
    print(f"Selected {len(selected_submissions)} submissions for processing.")
    with open(output_batch_file, "w") as f:
        if not selected_submissions:
            print("No new or continuable submissions found.")
        for submission in selected_submissions:
            origin = submission["origin"]
            submitter_id, inbox = origin["submitter_id"], origin["inbox"]
            submission_id = submission["submission_id"]
            target_path = f"results/{submitter_id}/{inbox}/{submission_id}/processed"
            f.write(f"{target_path}\n")
