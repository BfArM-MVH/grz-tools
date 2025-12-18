import json
import os
import sys
from contextlib import redirect_stderr, redirect_stdout

sys.path.append(os.path.dirname(__file__))
import shared

if snakemake.params.s3_access_key:
    os.environ["GRZ_S3__ACCESS_KEY"] = snakemake.params.s3_access_key
if snakemake.params.s3_secret:
    os.environ["GRZ_S3__SECRET"] = snakemake.params.s3_secret

inbox_config = snakemake.input.inbox_config_path
output_file = snakemake.output.submissions
submitter_id = snakemake.wildcards.submitter_id
inbox = snakemake.wildcards.inbox

stdout_log = snakemake.log.stdout
stderr_log = snakemake.log.stderr

with redirect_stdout(open(stdout_log, "w")), redirect_stderr(open(stderr_log, "w")):
    submissions_in_s3 = shared.scan_inbox_and_augment(inbox_config, submitter_id, inbox)
    with open(output_file, "w") as f:
        json.dump(submissions_in_s3, f, indent=2)
