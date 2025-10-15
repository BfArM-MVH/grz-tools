import json
import os
import subprocess
import tarfile
import time
from pathlib import Path

import pytest
import requests

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

SNAKEMAKE_BASE_CMD = PIXI_RUN_PREFIX + [
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
def docker_compose_file(project_root: Path) -> str:
    """Returns the path to the test-specific docker-compose.yml file."""
    return str(project_root / "packages/grz-watchdog/tests/docker-compose.test.yaml")


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
def test_environment(docker_compose_file: str):
    """Starts/stops container environment."""
    try:
        print("Starting TEST container environment in detached mode…")
        subprocess.run(
            [*CONTAINER_COMPOSE_CMD, "-f", docker_compose_file, "up", "--detach", "--build"],
            check=True,
            capture_output=True,
        )

        wait_for_services_ready(
            docker_compose_file, services_to_check=[GRZ_WATCHDOG_SERVICE_NAME, GRZ_SUBMITTER_SERVICE_NAME]
        )

        print("All services are up and ready for testing.")

        yield
    finally:
        print("\nStopping TEST container environment…")
        subprocess.run([*CONTAINER_COMPOSE_CMD, "-f", docker_compose_file, "down", "--volumes"], check=True)


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
