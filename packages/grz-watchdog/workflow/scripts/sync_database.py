import json
import os
import sys
from contextlib import redirect_stderr, redirect_stdout
from typing import Any

sys.path.append(os.path.dirname(__file__))
import shared
from shared import error_print, log_print


def register_submissions_with_db(submissions_json_list, db_config_path):  # noqa: C901
    try:
        db_states_info = shared.get_db_states(db_config_path)
    except Exception as e:
        error_print("Could not fetch current states from the database.")
        if hasattr(e, "stderr"):
            error_print(e.stderr)
        raise e

    current_db_states = {sub_id: info["state"] for sub_id, info in db_states_info.items()}

    all_s3_submissions = gather_submissions(submissions_json_list)

    if not all_s3_submissions:
        log_print("No submissions found in S3 input files. Nothing to do.")
        return []

    available_submissions = []
    for submission in all_s3_submissions:
        submission_id = submission["submission_id"]
        s3_state = submission["state"]

        if s3_state == "complete":
            target_db_state = "uploaded"
        elif s3_state == "incomplete":
            target_db_state = "uploading"
        else:
            log_print(f"Skipping submission {submission_id} with state {s3_state}.")
            continue

        current_db_state = current_db_states.get(submission_id)

        match (current_db_state, target_db_state):
            # A new, incomplete submission is found.
            case None, "uploading":
                log_print(f"Submission {submission_id} is new. Adding and setting state to 'uploading'.")
                shared.add_submission_to_db(db_config_path, submission_id)
                shared.update_submission_state_in_db(db_config_path, submission_id, "uploading")

            # A new, already complete submission is found.
            case None, "uploaded":
                log_print(f"Submission {submission_id} is new. Adding and setting state to 'uploaded'.")
                shared.add_submission_to_db(db_config_path, submission_id)
                shared.update_submission_state_in_db(db_config_path, submission_id, "uploaded")
                available_submissions.append(submission)

            # An existing incomplete submission is now complete.
            case "uploading", "uploaded":
                log_print(f"Updating state for {submission_id} from 'uploading' to 'uploaded'.")
                shared.update_submission_state_in_db(db_config_path, submission_id, "uploaded")
                available_submissions.append(submission)

            # "no-op" cases for clarity, nothing to be done here
            case "uploading", "uploading" | "uploaded", "uploaded":
                continue

            # catch-all for any other transition, which should be skipped.
            case from_state, to_state:
                log_print(
                    f"Skipping update for {submission_id}: No transition defined from '{from_state}' to '{to_state}'."
                )
                continue

    return available_submissions


def gather_submissions(submissions_json_list) -> list[Any]:
    all_s3_submissions = []

    for file_path in submissions_json_list:
        try:
            with open(file_path) as f:
                all_s3_submissions.extend(json.load(f))
        except (OSError, json.JSONDecodeError) as e:
            error_print(f"Error reading or parsing {file_path}.")
            raise e
    return all_s3_submissions


submissions_json_list = snakemake.input.submissions
db_config_path = snakemake.input.db_config_path
output_file = snakemake.output.submissions
stdout_log_path = snakemake.log.stdout
stderr_log_path = snakemake.log.stderr

with redirect_stdout(open(stdout_log_path, "w")), redirect_stderr(open(stderr_log_path, "w")):
    available_submissions = register_submissions_with_db(submissions_json_list, db_config_path)
    with open(output_file, "w") as f:
        json.dump(available_submissions, f, indent=2)
