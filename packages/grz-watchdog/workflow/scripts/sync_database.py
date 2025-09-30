import json
import subprocess
import sys
from contextlib import redirect_stderr, redirect_stdout
from os import PathLike

SUBPROCESS_TIMEOUT = 120


def log_print(*args, **kwargs):
    print(*args, **kwargs, file=sys.stdout)
    sys.stdout.flush()


def error_print(*args, **kwargs):
    print(*args, **kwargs, file=sys.stderr)
    sys.stderr.flush()


def get_current_db_states(db_config_path: PathLike | str) -> dict[str, str]:
    """
    Fetches the latest state for all submissions currently in the database.

    Returns:
        A dictionary mapping submission_id to its latest state (string).
    """
    log_print("Fetching current submission states from the database…")
    try:
        result = subprocess.run(
            ["grzctl", "db", "--config-file", db_config_path, "list", "--json"],
            check=True,
            text=True,
            capture_output=True,
            timeout=SUBPROCESS_TIMEOUT,
        )
        db_submissions = json.loads(result.stdout)

        db_states = {}
        for sub in db_submissions:
            if latest_state := sub.get("latest_state"):
                db_states[sub["id"]] = latest_state.get("state", "").casefold()

        return db_states

    except (subprocess.CalledProcessError, json.JSONDecodeError, subprocess.TimeoutExpired) as e:
        error_print("Could not fetch current states from the database.")
        if hasattr(e, "stderr"):
            error_print(e.stderr)
        raise e


def register_submissions_with_db(submissions_json_list, db_config_path):
    current_db_states = get_current_db_states(db_config_path)
    all_s3_submissions = []

    for file_path in submissions_json_list:
        try:
            with open(file_path) as f:
                all_s3_submissions.extend(json.load(f))
        except (OSError, json.JSONDecodeError) as e:
            error_print(f"Error reading or parsing {file_path}.")
            raise e

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

        current_state = current_db_states.get(submission_id)
        # A submission can be new (then its state is None) or in the process of being uploaded.
        allowed_previous_states = {None, "uploading"}

        if current_state not in allowed_previous_states:
            log_print(
                f"Skipping update for {submission_id}: "
                f"Its current DB state '{current_state}' is not in {allowed_previous_states}."
            )
            continue

        if current_state is None:
            try:
                log_print(f"Adding new submission {submission_id} to database…")
                subprocess.run(
                    ["grzctl", "db", "--config-file", db_config_path, "submission", "add", submission_id],
                    check=True,
                    text=True,
                    stdout=sys.stdout,
                    stderr=subprocess.PIPE,
                    timeout=SUBPROCESS_TIMEOUT,
                )
                log_print(f"Successfully added {submission_id}.")
            except subprocess.CalledProcessError as e:
                if "Duplicate submission ID" in e.stderr:
                    log_print(f"Submission {submission_id} already exists. Skipping add.")
                else:
                    error_print(e.stderr)
                    raise e
            except subprocess.TimeoutExpired as e:
                error_print(
                    f"Timeout: 'grzctl submission add' for {submission_id} exceeded {SUBPROCESS_TIMEOUT} seconds."
                )
                if e.stderr:
                    error_print(e.stderr)
                raise e

        try:
            log_print(f"Updating state for {submission_id} to '{target_db_state}'…")
            result = subprocess.run(
                [
                    "grzctl",
                    "db",
                    "--config-file",
                    db_config_path,
                    "submission",
                    "update",
                    submission_id,
                    target_db_state,
                ],
                check=True,
                text=True,
                capture_output=True,
                timeout=SUBPROCESS_TIMEOUT,
            )
            log_print(result.stdout)

            available_submissions.append(submission)
            log_print(f"Updated state for {submission_id} to {target_db_state}.")
        except subprocess.CalledProcessError as e:
            prompt_text = "Submission is currently in an 'Error' state"
            if prompt_text in e.stderr:
                error_print(
                    f"The database state for {submission_id} cannot be updated from 'Error' in non-interactive mode."
                )
                error_print(f"Captured error:\n{e.stderr}")
                raise e
            else:
                error_print(e.stderr)
                error_print(f"Error updating state for {submission_id}.")
                raise e
        except subprocess.TimeoutExpired as e:
            error_print(
                f"Timeout: 'grzctl submission update' for {submission_id} exceeded {SUBPROCESS_TIMEOUT} seconds."
            )
            if e.stderr:
                error_print(e.stderr)
            raise e

    return available_submissions


submissions_json_list = snakemake.input.submissions
db_config_path = snakemake.input.db_config_path
output_file = snakemake.output.submissions
stdout_log_path = snakemake.log.stdout
stderr_log_path = snakemake.log.stderr

with redirect_stdout(open(stdout_log_path, "w")), redirect_stderr(open(stderr_log_path, "w")):
    available_submissions = register_submissions_with_db(submissions_json_list, db_config_path)
    with open(output_file, "w") as f:
        json.dump(available_submissions, f, indent=2)
