import json
import os
import sys
from contextlib import redirect_stderr, redirect_stdout
from datetime import datetime

sys.path.append(os.path.dirname(__file__))
import shared


def select_submissions(submissions_path, db_config_path, limit):
    db_states = shared.get_db_states(db_config_path)

    with open(submissions_path) as f:
        synced_submissions = json.load(f)

    pending_submissions = []
    for submission in synced_submissions:
        submission_id = submission["submission_id"]
        latest_db_state = db_states.get(submission_id)

        if latest_db_state and (
            latest_db_state["state"] == "uploaded" or latest_db_state["state"] in shared.CONTINUABLE_STATES
        ):
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

with redirect_stdout(open(stdout_log_path, "w")), redirect_stderr(open(stderr_log_path, "w")):
    selected_submissions = select_submissions(submissions_path, db_config_path, limit)
    print(f"Selected {len(selected_submissions)} submissions for processing.")
    with open(output_batch_file, "w") as f:
        for submission in selected_submissions:
            origin = submission["origin"]
            submitter_id, inbox = origin["submitter_id"], origin["inbox"]
            submission_id = submission["submission_id"]
            target_path = f"results/{submitter_id}/{inbox}/{submission_id}/processed"
            f.write(f"{target_path}\n")
