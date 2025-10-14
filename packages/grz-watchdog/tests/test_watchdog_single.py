import re
import subprocess
from pathlib import Path

import pytest

from .conftest import (
    BUCKET_CONSENTED,
    BUCKET_INBOX,
    BUCKET_NONCONSENTED,
    CONTAINER_COMPOSE_CMD,
    CONTAINER_RUNTIME,
    GRZ_SUBMITTER_CONTAINER_NAME,
    GRZ_SUBMITTER_SERVICE_NAME,
    GRZ_WATCHDOG_CONTAINER_NAME,
    GRZ_WATCHDOG_SERVICE_NAME,
    INBOX,
    MINIO_SERVICE_NAME,
    PIXI_RUN_PREFIX,
    SNAKEMAKE_BASE_CMD,
    SUBMITTER_ID,
)

TEST_CASES = [
    ("panel", "123456789_2024-11-08_d0f805c5"),
    # ("wgs", "123456789_2024-10-28_e1bab61b"),
    ("wgs_lr", "123456789_2024-10-28_e1bab61b"),
]


def run_in_container(
    docker_compose_file: str, *args: str, service=GRZ_WATCHDOG_SERVICE_NAME
) -> subprocess.CompletedProcess:
    """Helper to execute a command inside the container via the compose service name."""
    command = [*CONTAINER_COMPOSE_CMD, "-f", docker_compose_file, "exec", "-T", service, *args]
    return subprocess.run(command, capture_output=True, text=True, check=True)


@pytest.fixture(scope="function")
def setup_and_submit(docker_compose_file: str, test_data_dir: Path, request):
    """A fixture to set up S3 and submit a specific test dataset."""
    submission_type = request.param
    local_data_path = test_data_dir / submission_type
    container_data_path = f"/tmp/{submission_type}"

    run_in_container(docker_compose_file, "sh", "-c", "rm -rf /workdir/results/* /tmp/*")
    run_in_container(docker_compose_file, "rm", "-f", "/workdir/results", service=GRZ_WATCHDOG_SERVICE_NAME)

    buckets_to_clean = [BUCKET_INBOX, BUCKET_CONSENTED, BUCKET_NONCONSENTED]
    for bucket in buckets_to_clean:
        cleanup_command = f"mc rm --recursive --force adm/{bucket}/* || true"
        run_in_container(docker_compose_file, "sh", "-c", cleanup_command, service=MINIO_SERVICE_NAME)

    run_in_container(
        docker_compose_file, "mc", "mb", "--ignore-existing", f"adm/{BUCKET_INBOX}", service=MINIO_SERVICE_NAME
    )

    subprocess.run(
        [CONTAINER_RUNTIME, "cp", str(local_data_path), f"{GRZ_SUBMITTER_CONTAINER_NAME}:{container_data_path}"],
        check=True,
    )

    run_in_container(
        docker_compose_file,
        *PIXI_RUN_PREFIX,
        "grz-cli",
        "submit",
        "--submission-dir",
        container_data_path,
        "--config-file",
        "/workdir/config/grz-cli.config.yaml",
        service=GRZ_SUBMITTER_SERVICE_NAME,
    )

    yield


@pytest.mark.parametrize("setup_and_submit, submission_id", TEST_CASES, indirect=["setup_and_submit"])
def test_single_valid_submission(docker_compose_file: str, setup_and_submit, submission_id: str):
    """Tests the successful end-to-end processing of a single valid submission."""
    print(f"\nRunning Snakemake for ({submission_id})...")
    target_file = f"results/{SUBMITTER_ID}/{INBOX}/{submission_id}/processed/without_qc"

    try:
        run_in_container(docker_compose_file, *SNAKEMAKE_BASE_CMD, "--cores", "1", target_file)
    except subprocess.CalledProcessError as e:
        print("\n--- SNAKEMAKE FAILED ---")
        print("Parsing stderr to find and dump specific log files...")

        log_path_pattern = re.compile(r"logs/[\w/.-]+\.log")
        found_logs = sorted(list(set(log_path_pattern.findall(e.stderr))))

        if not found_logs:
            print("\nCould not automatically find any log file paths in Snakemake's stderr.")
            print("--- Raw Snakemake stderr: ---")
            print(e.stderr)
        else:
            print(f"Found {len(found_logs)} log file(s) to dump: {', '.join(found_logs)}")
            for log_path in found_logs:
                try:
                    full_path_in_container = f"/workdir/{log_path}"
                    log_content = run_in_container(
                        docker_compose_file, "cat", full_path_in_container, service=GRZ_WATCHDOG_SERVICE_NAME
                    )
                    print(f"\n--- CONTENTS OF {log_path} ---")
                    print(log_content.stdout)
                    if log_content.stderr:
                        print(f"--- (stderr while cat'ing {log_path}) ---\n{log_content.stderr}")

                except subprocess.CalledProcessError as log_e:
                    print(f"\n--- Could not retrieve {log_path} ---")
                    print(log_e.stderr)

        pytest.fail(f"Snakemake workflow failed. See dumped log file contents above for details.", pytrace=False)

    print("Verifying outcomes...")
    db_result = run_in_container(
        docker_compose_file,
        *PIXI_RUN_PREFIX,
        "grzctl",
        "db",
        "--config-file",
        "/config/db.yaml",
        "submission",
        "show",
        submission_id,
    )
    assert "State: Finished" in db_result.stdout, "Submission state was not 'Finished' in the database."

    inbox_path = f"adm/{BUCKET_INBOX}/{submission_id}"
    inbox_ls_result = run_in_container(
        docker_compose_file, "mc", "ls", "--recursive", inbox_path, service=MINIO_SERVICE_NAME
    )
    inbox_files = {line.split()[-1] for line in inbox_ls_result.stdout.strip().split("\n")}
    assert inbox_files == {f"{inbox_path}/cleaned", f"{inbox_path}/metadata/metadata.json"}, (
        "Inbox was not cleaned correctly."
    )

    archive_path = f"adm/{BUCKET_CONSENTED}/{submission_id}"
    archive_ls_result = run_in_container(docker_compose_file, "mc", "ls", archive_path, service=MINIO_SERVICE_NAME)
    assert "metadata/" in archive_ls_result.stdout, "Archived submission missing metadata directory."
    assert "files/" in archive_ls_result.stdout, "Archived submission missing files directory."
