import hashlib
import json
import os
import re
import shutil
import subprocess
import tarfile
import time
from pathlib import Path

import pytest
import requests

from grz_pydantic_models.submission.metadata import GrzSubmissionMetadata

PROJECT_ROOT = Path(__file__).parent.parent.parent.parent
DOCKER_COMPOSE_FILE = str(PROJECT_ROOT / "packages/grz-watchdog/tests/docker-compose.test.yaml")

CONTAINER_RUNTIME = "podman"  # or "docker"
CONTAINER_COMPOSE_CMD = ["podman-compose"]  # or ["docker", "compose"]
GRZ_WATCHDOG_CONTAINER_NAME = "grz-watchdog-test-runner"
GRZ_WATCHDOG_SERVICE_NAME = "grz-watchdog"
GRZ_SUBMITTER_SERVICE_NAME = "grz-submitter"
GRZ_SUBMITTER_CONTAINER_NAME = "grz-submitter-test-runner"
MINIO_SERVICE_NAME = "minio"

SUBMITTER_ID = "123456789"
INBOX = "test1"
BUCKET_INBOX = "test1"
BUCKET_CONSENTED = "consented"
BUCKET_NONCONSENTED = "nonconsented"

PIXI_RUN_PREFIX = [
    "pixi",
    "run",
    "--manifest-path",
    "/workspace/packages/grz-watchdog/tests/pixi.toml",
    "--",
]

SNAKEMAKE_BASE_CMD = [
    *PIXI_RUN_PREFIX,
    "snakemake",
    "--workflow-profile",
    "/workspace/packages/grz-watchdog/workflow/profiles/test",
    "--snakefile",
    "/workspace/packages/grz-watchdog/workflow/Snakefile",
    "--directory",
    "/workdir",
    "--conda-prefix",
    "/conda-envs",
]


@pytest.fixture(scope="session")
def project_root() -> Path:
    """Returns the grz-tools root directory."""
    return Path(__file__).parent.parent.parent.parent


@pytest.fixture(scope="session")
def test_data_dir(tmpdir_factory, version: str = "0.2.2"):
    """
    Downloads (with caching) and extracts all test data into a single session directory.
    Each dataset is placed in its own named subdirectory.
    """
    cache_dir = Path(__file__).parent / ".test_cache" / version
    cache_dir.mkdir(parents=True, exist_ok=True)

    session_data_dir = tmpdir_factory.mktemp("test_data_session")

    for what in ("panel", "references", "wgs", "wgs_lr"):
        cached_tarball_path = cache_dir / f"{what}.tgz"

        if not cached_tarball_path.exists():
            print(f"\nDownloading test data to {cached_tarball_path}…")
            url = f"https://github.com/twrightsman/grz-mini-test-data/releases/download/version/{version}/{what}.tgz"
            with requests.get(url, stream=True) as r:
                r.raise_for_status()
                with open(cached_tarball_path, "wb") as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)
        else:
            print(f"\nUsing cached test data from {cached_tarball_path}…")

        extract_dir = session_data_dir / what
        with tarfile.open(cached_tarball_path, "r:gz") as tar:
            tar.extractall(path=extract_dir, filter="data")

    return session_data_dir


# workaround for missing "--wait" in podman: https://github.com/containers/podman-compose/issues/710
def wait_for_services_ready(docker_compose_file: str, services_to_check: list[str], timeout: int = 120):
    """
    Polls `podman-compose ps` until all specified services are in a 'running' state.
    """
    print(f"Waiting for services to become fully running: {services_to_check}")
    start_time = time.time()
    not_ready_yet = set(services_to_check)

    while time.time() - start_time < timeout:
        try:
            ps_result = subprocess.run(
                [*CONTAINER_COMPOSE_CMD, "-f", docker_compose_file, "ps", "--format", "json"],
                check=True,
                capture_output=True,
                text=True,
            )
            services_status = json.loads(ps_result.stdout)

            for service_state in services_status:
                service_name = service_state.get("Labels", {}).get("com.docker.compose.service")

                if service_name in not_ready_yet and service_state.get("State") == "running":
                    print(f"Service '{service_name}' is running.")
                    not_ready_yet.remove(service_name)

            if not not_ready_yet:
                print("All specified services are up and running.")
                return

            time.sleep(2)
        except (subprocess.CalledProcessError, json.JSONDecodeError, KeyError) as e:
            print(f"Waiting for services… (Encountered temporary error: {e})")
            time.sleep(2)

    pytest.fail(f"The following services did not become ready within the {timeout}s timeout: {list(not_ready_yet)}")


@pytest.fixture(scope="session", autouse=True)
def test_environment():
    """Starts/stops container environment."""
    try:
        print("Starting TEST container environment in detached mode…")
        subprocess.run(
            [*CONTAINER_COMPOSE_CMD, "-f", DOCKER_COMPOSE_FILE, "up", "--detach", "--build"],
            check=True,
            capture_output=True,
        )

        wait_for_services_ready(
            DOCKER_COMPOSE_FILE, services_to_check=[GRZ_WATCHDOG_SERVICE_NAME, GRZ_SUBMITTER_SERVICE_NAME]
        )

        print("All services are up and ready for testing.")

        yield
    finally:
        print("Stopping TEST container environment…")
        subprocess.run([*CONTAINER_COMPOSE_CMD, "-f", DOCKER_COMPOSE_FILE, "down", "--volumes"], check=True)


@pytest.fixture(scope="session", autouse=True)
def mock_api_certificates(project_root: Path):
    """
    Generates self-signed certificates for the mock API service if they don't exist.
    """
    cert_dir = project_root / "packages/grz-watchdog/tests/generated"
    key_file = cert_dir / "key.pem"
    cert_file = cert_dir / "cert.pem"

    if key_file.exists() and cert_file.exists():
        print("Mock API certificates already exist. Skipping generation.")
        return

    print(f"Generating Mock API certificates in {cert_dir}…")
    os.makedirs(cert_dir, exist_ok=True)

    try:
        subprocess.run(
            [
                "openssl",
                "req",
                "-x509",
                "-newkey",
                "rsa:4096",
                "-keyout",
                str(key_file),
                "-out",
                str(cert_file),
                "-sha256",
                "-days",
                "365",
                "-nodes",
                "-subj",
                "/CN=bfarm-mock-api",
                "-addext",
                "subjectAltName = DNS:bfarm-mock-api",
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        print("Certificates generated successfully.")
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        pytest.fail(f"Failed to generate self-signed certificates with openssl.\nError: {e.stderr}", pytrace=False)


def run_in_container(*args: str, service: str = GRZ_WATCHDOG_SERVICE_NAME) -> subprocess.CompletedProcess:
    """Helper to execute a command inside the container via the compose service name."""
    command = [*CONTAINER_COMPOSE_CMD, "-f", DOCKER_COMPOSE_FILE, "exec", "-T", service, *args]
    return subprocess.run(command, capture_output=True, text=True, check=True)


class BaseTest:
    """A base class for tests, containing reusable helper methods."""

    def _submit_data(self, local_data_path: Path):
        """Helper to submit a specific test dataset."""
        local_data_path = Path(local_data_path)
        container_data_path = f"/tmp/{local_data_path.name}"
        print(f"\n--- Submitting Data from {local_data_path.name} ---")

        subprocess.run(
            [CONTAINER_RUNTIME, "cp", str(local_data_path), f"{GRZ_SUBMITTER_CONTAINER_NAME}:{container_data_path}"],
            check=True,
        )
        run_in_container(
            *PIXI_RUN_PREFIX,
            "grz-cli",
            "submit",
            "--submission-dir",
            container_data_path,
            "--config-file",
            "/workdir/config/grz-cli.config.yaml",
            service=GRZ_SUBMITTER_SERVICE_NAME,
        )

    def _handle_watchdog_failure(self, e: subprocess.CalledProcessError):
        """Dump logs from a failed grz-watchdog run."""
        log_path_pattern = re.compile(r"logs/[\w/.-]+\.log")
        found_logs = sorted(list(set(log_path_pattern.findall(e.stderr))))
        print(e.stderr)

        if not found_logs:
            print("\nCould not automatically find any log file paths in Snakemake's stderr.")
            return

        print(f"Found {len(found_logs)} log file(s) to dump: {', '.join(found_logs)}")
        for log_path in found_logs:
            try:
                full_path = f"/workdir/{log_path}"
                log_content = run_in_container("cat", full_path)
                print(f"\n--- CONTENTS OF {log_path} ---")
                print(log_content.stdout)
            except subprocess.CalledProcessError as log_e:
                print(f"\n--- Could not retrieve {log_path} ---\n{log_e.stderr}")

    def _run_watchdog(self, target: str, cores: int = 1):
        """Run grz-watchdog and handle failures."""
        print(f"\nRunning grz-watchdog for target: {target}…")
        try:
            run_in_container(*SNAKEMAKE_BASE_CMD, "--cores", str(cores), target)
        except subprocess.CalledProcessError as e:
            self._handle_watchdog_failure(e)
            pytest.fail("grz-watchdog failed. See dumped log contents above for details.", pytrace=False)

    def _run_watchdog_expect_fail(self, target: str, cores: int = 1):
        """Run grz-watchdog and assert that it fails."""
        print(f"\nRunning grz-watchdog for target: {target} (expecting failure)…")
        with pytest.raises(subprocess.CalledProcessError) as excinfo:
            run_in_container(*SNAKEMAKE_BASE_CMD, "--cores", str(cores), target)

        self._handle_watchdog_failure(excinfo.value)
        return excinfo.value

    def start_background_process(self, command: list[str]) -> subprocess.Popen:
        """Starts a background process inside the container."""
        full_command = [*CONTAINER_COMPOSE_CMD, "-f", DOCKER_COMPOSE_FILE, "exec", GRZ_WATCHDOG_SERVICE_NAME, *command]
        process = subprocess.Popen(
            full_command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            universal_newlines=True,
        )
        return process

    def stop_background_process(self, process: subprocess.Popen):
        """Stops a running background process."""
        print("Stopping background process...")
        try:
            # Send TERM signal for graceful shutdown
            run_in_container("pkill", "-SIGTERM", "-f", "snakemake")
            stdout, _ = process.communicate(timeout=120)
            print(stdout)
        except (subprocess.TimeoutExpired, subprocess.CalledProcessError):
            print("Process did not terminate gracefully, killing.")
            try:
                run_in_container("pkill", "-SIGKILL", "-f", "snakemake")
            except subprocess.CalledProcessError:
                pass
            process.kill()
            stdout, _ = process.communicate()
            print(stdout)

    def _verify_db_state(self, submission_id: str, expected_state: str):
        """Verifies the final state of a submission in the database."""
        print(f"Verifying database state for {submission_id}…")
        result = run_in_container(
            *PIXI_RUN_PREFIX, "grzctl", "db", "--config-file", "/config/configs/db.yaml", "list", "--json"
        )
        submissions = json.loads(result.stdout)
        target = next((s for s in submissions if s.get("id") == submission_id), None)

        assert target is not None, f"Submission '{submission_id}' not found in the database."
        final_state = target.get("latest_state", {}).get("state")
        assert final_state == expected_state, (
            f"Submission state was not '{expected_state}'. Actual state: '{final_state}'"
        )
        print(f"OK: Database state is '{final_state}'.")

    def _verify_inbox_cleaned(self, submission_id: str):
        """Verifies that the S3 inbox was cleaned correctly."""
        print(f"Verifying inbox state for {submission_id}…")

        _ = run_in_container(
            *PIXI_RUN_PREFIX,
            "mc alias set verify http://minio:9000 minioadmin minioadmin",
            service=GRZ_WATCHDOG_SERVICE_NAME,
        )
        result = run_in_container(
            *PIXI_RUN_PREFIX,
            f"mc ls --recursive verify/{BUCKET_INBOX}/{submission_id}",
            service=GRZ_WATCHDOG_SERVICE_NAME,
        )

        files = {line.split()[-1] for line in result.stdout.strip().split("\n")}
        expected_files = {"cleaned", "metadata/metadata.json"}

        assert files == expected_files, f"Inbox was not cleaned correctly. Found: {files}"
        print("OK: Inbox is cleaned correctly.")

    def _verify_archived(self, submission_id: str, bucket: str):
        """Verifies that a submission was correctly archived to the target bucket."""
        print(f"Verifying archive state for {submission_id} in bucket {bucket}...")
        archive_path = f"verify/{bucket}/{submission_id}"

        _ = run_in_container(
            *PIXI_RUN_PREFIX,
            "mc alias set verify http://minio:9000 minioadmin minioadmin",
            service=GRZ_WATCHDOG_SERVICE_NAME,
        )
        result = run_in_container(*PIXI_RUN_PREFIX, f"mc ls {archive_path}", service=GRZ_WATCHDOG_SERVICE_NAME)

        assert "metadata/" in result.stdout, "Archived submission missing metadata directory: " + result.stdout
        assert "files/" in result.stdout, "Archived submission missing files directory: " + result.stdout
        print(f"OK: Submission correctly archived in {bucket}.")


@pytest.fixture(scope="class")
def setup_class_environment():
    """A fixture to clean the environment once before any tests in the class run."""
    print("\n--- Cleaning S3 buckets and workspace before tests ---")
    run_in_container("sh", "-c", "rm -rf /workdir/results/* /tmp/*")

    buckets_to_clean = [BUCKET_INBOX, BUCKET_CONSENTED, BUCKET_NONCONSENTED]
    for bucket in buckets_to_clean:
        cleanup_command = f"mc rm --recursive --force adm/{bucket}/* || true"
        run_in_container("sh", "-c", cleanup_command, service=MINIO_SERVICE_NAME)

    run_in_container("mc", "mb", "--ignore-existing", f"adm/{BUCKET_INBOX}", service=MINIO_SERVICE_NAME)
    print("\n--- Finished cleaning ---")
    yield


def _create_variant_submission(base_dir: Path, variant_name_suffix: str, tmp_path: Path):
    """Create a modified submission with a different tanG."""
    base_dir, tmp_dir = Path(base_dir), Path(tmp_path)
    variant_dir = tmp_path / f"{base_dir.name}_{variant_name_suffix}"
    if variant_dir.exists():
        shutil.rmtree(variant_dir)
    shutil.copytree(base_dir, variant_dir)

    metadata_path = variant_dir / "metadata" / "metadata.json"
    with open(metadata_path) as f:
        metadata = GrzSubmissionMetadata(**json.load(f))

    original_tang = metadata.submission.tan_g
    new_tang = hashlib.sha256(f"{original_tang}{variant_name_suffix}".encode()).hexdigest()
    metadata.submission.tan_g = new_tang

    with open(metadata_path, "w") as f:
        f.write(metadata.model_dump_json(by_alias=True))

    submission_id = metadata.submission_id
    print(metadata)

    return variant_dir, submission_id
