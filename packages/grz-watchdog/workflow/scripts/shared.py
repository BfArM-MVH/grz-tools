import json
import os
import subprocess
import sys

sys.path.append(os.path.dirname(__file__))

SUBPROCESS_TIMEOUT = 120

CONTINUABLE_STATES = {
    "downloaded",
    "decrypted",
    "validated",
    "encrypted",
    "archived",
    "reported",
    "qced",
    "cleaned",
}


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


def scan_inbox_and_augment(grzctl_config, submitter_id, inbox):
    """
    Scans a single inbox using 'grzctl list' and augments submission data with its origin (submitter_id, inbox).
    """
    result = run_grzctl_command(
        [
            "--config",
            grzctl_config,
            "list",
            "--submitter-id",
            submitter_id,
            "--inbox",
            inbox,
            "--json",
            "--show-cleaned",
            "--limit",
            "1000000",
        ]
    )
    submissions = json.loads(result.stdout)
    for submission in submissions:
        submission["origin"] = {"submitter_id": submitter_id, "inbox": inbox}
    return submissions


def get_db_states(grzctl_config):
    """
    Fetches all submissions from the db and returns their latest states and timestamps.
    """
    result = run_grzctl_command(["--config", grzctl_config, "db", "list", "--json", "--limit", "1000000"])
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
