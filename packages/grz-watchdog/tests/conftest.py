import subprocess
import tarfile
import time
from pathlib import Path

import pytest
import requests

CONTAINER_RUNTIME = "podman"  # or "docker"
CONTAINER_COMPOSE_CMD = ["podman-compose"]  # or ["docker", "compose"]
CONTAINER_NAME = "grz-watchdog-test-runner"
SERVICE_NAME = "grz-watchdog"

SUBMITTER_ID = "123456789"
INBOX = "test1"
BUCKET_INBOX = "adm/test1"
BUCKET_CONSENTED = "adm/consented"
BUCKET_NONCONSENTED = "adm/nonconsented"


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
            print(f"\nDownloading test data to {cached_tarball_path}...")
            url = f"https://github.com/twrightsman/grz-mini-test-data/releases/download/version/{version}/{what}.tgz"
            with requests.get(url, stream=True) as r:
                r.raise_for_status()
                with open(cached_tarball_path, "wb") as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)
        else:
            print(f"\nUsing cached test data from {cached_tarball_path}...")

        extract_path = session_data_dir / what
        extract_path.mkdir()
        with tarfile.open(cached_tarball_path, "r:gz") as tar:
            tar.extractall(path=extract_path)

    return session_data_dir


@pytest.fixture(scope="session", autouse=True)
def test_environment(docker_compose_file: str):
    """Starts/stops container environment and waits for it to be ready."""
    try:
        print("Starting TEST container environment in detached mode...")

        subprocess.Popen(
            [*CONTAINER_COMPOSE_CMD, "-f", docker_compose_file, "up", "--detach", "--build"],
        )

        print("Waiting for services to become healthy and configured...")
        timeout = 300
        start_time = time.time()
        mc_alias_set = False
        while time.time() - start_time < timeout:
            try:
                subprocess.run(
                    [
                        *CONTAINER_COMPOSE_CMD,
                        "-f",
                        docker_compose_file,
                        "exec",
                        "-T",
                        SERVICE_NAME,
                        "mc",
                        "alias",
                        "set",
                        "adm",
                        "http://minio:9000",
                        "minioadmin",
                        "minioadmin",
                    ],
                    check=True,
                    capture_output=True,
                )
                print("Services are ready and mc alias is configured.")
                mc_alias_set = True
                break
            except subprocess.CalledProcessError as e:
                print("Services not ready yet, retrying in 10 seconds...", e)
                time.sleep(10)

        if not mc_alias_set:
            pytest.fail(f"Services did not become ready within the {timeout} second timeout.")

        yield
    finally:
        print("\nStopping TEST container environment...")
        subprocess.run([*CONTAINER_COMPOSE_CMD, "-f", docker_compose_file, "down", "--volumes"], check=True)
