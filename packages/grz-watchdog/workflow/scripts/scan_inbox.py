import json
import subprocess
import os
from contextlib import redirect_stdout, redirect_stderr

if snakemake.params.s3_access_key:
    os.environ["GRZ_S3__ACCESS_KEY"] = snakemake.params.s3_access_key
if snakemake.params.s3_secret:
    os.environ["GRZ_S3__SECRET"] = snakemake.params.s3_secret


def scan_inbox_and_augment(inbox_config, submitter_id, inbox):
    try:
        result = subprocess.run(
            ["grzctl", "list", "--config-file", inbox_config, "--json", "--show-cleaned"],
            capture_output=True,
            text=True,
            check=True,
        )
        submissions = json.loads(result.stdout)

        # explicitly include submitter_id and inbox in output entry
        for submission in submissions:
            submission["origin"] = {"submitter_id": submitter_id, "inbox": inbox}

        return submissions

    except (subprocess.CalledProcessError, json.JSONDecodeError) as e:
        raise e


inbox_config = snakemake.input.inbox_config_path
output_file = snakemake.output.submissions
submitter_id = snakemake.wildcards.submitter_id
inbox = snakemake.wildcards.inbox

stdout_log = snakemake.log.stdout
stderr_log = snakemake.log.stderr

with redirect_stdout(open(stdout_log, "w")), redirect_stderr(open(stderr_log, "w")):
    submissions_in_s3 = scan_inbox_and_augment(inbox_config, submitter_id, inbox)
    with open(output_file, "w") as f:
        json.dump(submissions_in_s3, f, indent=2)
