import json
import subprocess
import sys
from contextlib import redirect_stderr, redirect_stdout

SUBPROCESS_TIMEOUT = 120


def log_print(*args, **kwargs):
    print(*args, **kwargs, file=sys.stdout)
    sys.stdout.flush()


def error_print(*args, **kwargs):
    print(*args, **kwargs, file=sys.stderr)
    sys.stderr.flush()


def register_submissions_with_db(submissions_json_list, db_config_path):
    all_submissions = []

    for file_path in submissions_json_list:
        try:
            with open(file_path) as f:
                all_submissions.extend(json.load(f))
        except (OSError, json.JSONDecodeError) as e:
            raise e

    if not all_submissions:
        log_print("No submissions found. Nothing to do.")
        return []

    available_submissions = []
    for submission in all_submissions:
        submission_id = submission["submission_id"]
        s3_state = submission["state"]

        if s3_state == "complete":
            db_state = "uploaded"
        elif s3_state == "incomplete":
            db_state = "uploading"
        else:
            log_print(f"Skipping submission {submission_id} with state {s3_state}.")
            continue

        try:
            log_print(f"Adding submission {submission_id} to database…")
            subprocess.run(
                ["grzctl", "db", "--config-file", db_config_path, "submission", "add", submission_id],
                check=True,
                text=True,
                stdout=sys.stdout,
                stderr=subprocess.PIPE,
                timeout=SUBPROCESS_TIMEOUT,
            )
            log_print(f"Added submission {submission_id} to database.")
        except subprocess.CalledProcessError as e:
            if "Duplicate submission ID " in e.stderr:
                log_print(f"Submission {submission_id} already exists in database. Skipping.")
            else:
                error_print(e.stderr)
                raise e
        except subprocess.TimeoutExpired as e:
            error_print(f"Timeout: 'grzctl submission add' for {submission_id} exceeded {SUBPROCESS_TIMEOUT} seconds.")
            if e.stderr:
                error_print(e.stderr)
            raise e

        try:
            log_print(f"Updating state for {submission_id} to {db_state}…")
            subprocess.run(
                ["grzctl", "db", "--config-file", db_config_path, "submission", "update", submission_id, db_state],
                check=True,
                text=True,
                stdout=sys.stdout,
                stderr=subprocess.PIPE,
                timeout=SUBPROCESS_TIMEOUT,
            )
            available_submissions.append(submission)
            log_print(f"Updated state for {submission_id} to {db_state}.")
        except subprocess.CalledProcessError as e:
            error_print(e.stderr)
            error_print(f"Error updating state for {submission_id}.")
            raise e
        except subprocess.TimeoutExpired as e:
            error_print(
                f"Timeout: 'grzctl submission update' for {submission_id} exceeded {SUBPROCESS_TIMEOUT} seconds."
            )
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
