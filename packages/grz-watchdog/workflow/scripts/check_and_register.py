import os
import subprocess
import sys

sys.path.append(os.path.dirname(__file__))
import shared


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
db_config = snakemake.input.db_config_path
log_stdout = snakemake.log.stdout
log_stderr = snakemake.log.stderr

with open(log_stdout, "w") as f_out, open(log_stderr, "w") as f_err:
    print(f"Checking S3 status for submission '{submission_id}'...", file=f_out)
    submissions = shared.scan_inbox_and_augment(inbox_config=inbox_config, submitter_id=submitter_id, inbox=inbox)
    try:
        submission = next(filter(lambda s: s["submission_id"] == submission_id, submissions))
    except StopIteration:
        raise RuntimeError(f"Submission '{submission_id}' not found in S3 inbox.") from None

    s3_state = submission.get("state")
    print(f"Found S3 state: '{s3_state}'", file=f_out)

    if s3_state != "complete":
        raise RuntimeError(f"Submission '{submission_id}' is not ready for processing. S3 status is '{s3_state}'.")

    try:
        print(f"Registering submission '{submission_id}' in the database...", file=f_out)
        add_cmd = ["grzctl", "db", "--config-file", db_config, "submission", "add", submission_id]
        subprocess.run(add_cmd, capture_output=True, text=True)  # noqa: S603

        print("Setting database state to 'uploaded'...", file=f_out)
        update_cmd = ["grzctl", "db", "--config-file", db_config, "submission", "update", submission_id, "uploaded"]
        run_command(update_cmd)

        print("Submission successfully checked and registered.", file=f_out)

    except Exception as e:
        print(str(e), file=f_err)
        sys.exit(1)
