import json
import os
import subprocess
import sys

sys.path.append(os.path.dirname(__file__))

SUBPROCESS_TIMEOUT = 120


def log_print(*args, **kwargs):
    print(*args, **kwargs, file=sys.stdout)
    sys.stdout.flush()


def error_print(*args, **kwargs):
    print(*args, **kwargs, file=sys.stderr)
    sys.stderr.flush()


def run_grzctl_command(cmd, check=True):
    """Helper to run a grzctl command."""
    try:
        result = subprocess.run(  # noqa: S603
            ["grzctl", *cmd],  # noqa: S607
            check=check,
            text=True,
            capture_output=True,
            timeout=SUBPROCESS_TIMEOUT,
        )
        return result
    except (subprocess.CalledProcessError, json.JSONDecodeError, subprocess.TimeoutExpired) as e:
        raise e


def scan_inbox_and_augment(inbox_config, submitter_id, inbox):
    """
    Scans a single inbox using 'grzctl list' and augments submission data with its origin (submitter_id, inbox).
    """
    result = run_grzctl_command(
        ["list", "--config-file", inbox_config, "--json", "--show-cleaned", "--limit", "1000000"]
    )
    submissions = json.loads(result.stdout)
    for submission in submissions:
        submission["origin"] = {"submitter_id": submitter_id, "inbox": inbox}
    return submissions


def get_db_states(db_config_path):
    """
    Fetches all submissions from the db and returns their latest states and timestamps.
    """
    result = run_grzctl_command(["db", "--config-file", db_config_path, "list", "--json", "--limit", "1000000"])
    if not result.stdout:
        return {}

    db_states = {}
    for entry in json.loads(result.stdout):
        if latest_state := entry.get("latest_state"):
            db_states[entry["id"]] = {
                "state": latest_state.get("state", "").casefold(),
                "timestamp": latest_state.get("timestamp"),
            }
    return db_states


def add_submission_to_db(db_config_path, submission_id):
    """Runs 'grzctl db submission add {submission_id}'."""
    try:
        run_grzctl_command(["db", "--config-file", db_config_path, "submission", "add", submission_id])
    except subprocess.CalledProcessError as e:
        error_print(f"An unexpected error occurred for {submission_id}:")
        error_print(e.stderr)
        raise e
    except subprocess.TimeoutExpired as e:
        error_print(f"Timeout occurred while processing {submission_id}.")
        if e.stderr:
            error_print(e.stderr)


def update_submission_state_in_db(db_config_path, submission_id, state):
    """Runs 'grzctl db submission update {submission_id} {state}'."""
    try:
        run_grzctl_command(["db", "--config-file", db_config_path, "submission", "update", submission_id, state])
    except subprocess.CalledProcessError as e:
        if "Submission is currently in an 'Error' state" in e.stderr:
            error_print(f"The state for {submission_id} cannot be updated from 'Error' in non-interactive mode.")
            error_print(f"Captured error:\n{e.stderr}")
        else:
            error_print(f"An unexpected error occurred for {submission_id}:")
            error_print(e.stderr)
        raise e
    except subprocess.TimeoutExpired as e:
        error_print(f"Timeout occurred while updating {submission_id}.")
        if e.stderr:
            error_print(e.stderr)
        raise e
