import os
import subprocess
import sys

sys.path.append(os.path.dirname(__file__))
import shared
from shared import log_print


def run_command(cmd):
    try:
        return subprocess.run(cmd, check=True, text=True, capture_output=True, timeout=120)  # noqa: S603
    except subprocess.CalledProcessError as e:
        print(f"Command failed: {' '.join(cmd)}\nSTDERR:\n{e.stderr}", file=sys.stderr)
        raise


submission_id = snakemake.wildcards.submission_id
submitter_id = snakemake.wildcards.submitter_id
inbox = snakemake.wildcards.inbox
inbox_config = snakemake.input.inbox_config_path
db_config_path = snakemake.input.db_config_path
log_stdout = snakemake.log.stdout
log_stderr = snakemake.log.stderr

with open(log_stdout, "w") as f_out, open(log_stderr, "w") as f_err:
    print(f"Checking S3 status for submission '{submission_id}'...", file=f_out)
    submissions = shared.scan_inbox_and_augment(inbox_config=inbox_config, submitter_id=submitter_id, inbox=inbox)
    try:
        submission = next(filter(lambda s: s["submission_id"] == submission_id, submissions))
    except StopIteration:
        raise RuntimeError(f"Submission '{submission_id}' not found in S3 inbox.") from None

    db_state: str = (submission.get("database_state", "missing") or "missing").casefold()

    submission_id = submission["submission_id"]
    s3_state = submission["state"]

    if s3_state == "complete":
        target_db_state = "uploaded"
    elif s3_state == "incomplete":
        target_db_state = "uploading"
    else:
        log_print(f"Skipping submission {submission_id} with state {s3_state}.")
        sys.exit(0)

    current_db_state = None if db_state == "missing" else db_state

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

        # An existing incomplete submission is now complete.
        case "uploading", "uploaded":
            log_print(f"Updating state for {submission_id} from 'uploading' to 'uploaded'.")
            shared.update_submission_state_in_db(db_config_path, submission_id, "uploaded")

        # If we already set the state to "uploaded" in a previous run which subsequently got interrupted,
        # the submission is still available for processing
        case "uploaded", "uploaded":
            pass

        # "no-op" cases for clarity, nothing to be done here
        case "uploading", "uploading":
            pass

        # catch-all for any other transition, which should be skipped.
        case from_state, to_state:
            log_print(
                f"Skipping update for {submission_id}: No transition defined from '{from_state}' to '{to_state}'."
            )
            pass
