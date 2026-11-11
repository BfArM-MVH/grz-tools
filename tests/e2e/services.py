"""Services that grz-cli and grzctl talk to in the end-to-end tests."""

import datetime
import ipaddress
import json
import os
import shutil
import ssl
import subprocess
import threading
import time
from collections.abc import Iterator
from contextlib import contextmanager
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

import boto3
import botocore.exceptions
from cryptography import x509
from cryptography.hazmat.primitives import hashes, serialization
from cryptography.hazmat.primitives.asymmetric import ec
from cryptography.x509.oid import NameOID

SEAWEEDFS_IMAGE = "docker.io/chrislusf/seaweedfs:4.53"

ADMIN = ("e2e-admin", "e2e-admin-secret")
"""S3 credentials of the tests themselves, for setup and inspection."""

INBOX = "grz-inbox"
ARCHIVE_CONSENTED = "grz-archive-consented"
ARCHIVE_NON_CONSENTED = "grz-archive-non-consented"
ARCHIVE_INTERROGATION = "grz-archive-interrogation"
BUCKETS = (INBOX, ARCHIVE_CONSENTED, ARCHIVE_NON_CONSENTED, ARCHIVE_INTERROGATION)

_ARCHIVE = ("s3:GetObject", "s3:PutObject", "s3:AbortMultipartUpload", "s3:ListBucket")
PERMISSIONS: dict[str, dict[str, tuple[str, ...]]] = {
    "le": {INBOX: ("s3:GetObject", "s3:PutObject", "s3:AbortMultipartUpload", "s3:ListBucket")},
    "inbox": {INBOX: ("s3:GetObject", "s3:PutObject", "s3:DeleteObject", "s3:ListBucket")},
    "interrogation": {
        ARCHIVE_INTERROGATION: (
            "s3:GetObject",
            "s3:PutObject",
            "s3:DeleteObject",
            "s3:AbortMultipartUpload",
            "s3:ListBucket",
        )
    },
    "consented": {ARCHIVE_CONSENTED: _ARCHIVE, ARCHIVE_INTERROGATION: ("s3:GetObject", "s3:ListBucket")},
    "non_consented": {ARCHIVE_NON_CONSENTED: _ARCHIVE, ARCHIVE_INTERROGATION: ("s3:GetObject", "s3:ListBucket")},
}
"""S3 actions of each set of credentials, by bucket: the LE's in grz-cli, and those of the grzctl config.

The README of grzctl lists the same, and ``test_permissions_doc`` keeps both in line.
"""


def credentials(role: str) -> tuple[str, str]:
    """Return access key and secret of a role of :data:`PERMISSIONS`."""
    return f"e2e-{role}", f"e2e-{role}-secret"


def container_engine() -> str:
    """Return the container engine: ``GRZ_E2E_CONTAINER_ENGINE``, else docker or podman, whichever is installed."""
    if engine := os.environ.get("GRZ_E2E_CONTAINER_ENGINE"):
        return engine
    for engine in ("docker", "podman"):
        if shutil.which(engine):
            return engine
    raise RuntimeError("The end-to-end tests need docker or podman, or GRZ_E2E_CONTAINER_ENGINE.")


@contextmanager
def seaweedfs(work_dir: Path, engine_env: dict[str, str]) -> Iterator[str]:
    """Run SeaweedFS with an S3 gateway in a container, with the GRZ's buckets and a user per role of :data:`PERMISSIONS`.

    :param work_dir: Directory for the identities file.
    :param engine_env: Environment of the container engine; rootless podman needs the real ``HOME``.
    :returns: The S3 endpoint URL, on a free port of localhost.
    """
    identities = work_dir / "s3.json"
    users = [{"name": "admin", "credentials": [_keys(ADMIN)], "actions": ["Admin", "Read", "Write", "List", "Tagging"]}]
    # the roles get no actions of their own; the bucket policies grant them
    users += [
        {"name": credentials(role)[0], "credentials": [_keys(credentials(role))], "actions": []} for role in PERMISSIONS
    ]
    identities.write_text(json.dumps({"identities": users}))
    engine = container_engine()
    run = [engine, "run", "--rm", "--detach", "--publish", "127.0.0.1::8333"]
    run += ["--volume", f"{identities}:/etc/seaweedfs/s3.json:ro,z", SEAWEEDFS_IMAGE]
    run += ["server", "-dir=/data", "-s3", "-s3.config=/etc/seaweedfs/s3.json"]
    started = subprocess.run(run, env=engine_env, check=False, capture_output=True, text=True)
    if started.returncode:
        raise RuntimeError(f"{engine} run exited with {started.returncode}: {started.stderr}")
    container = started.stdout.strip()
    try:
        port = subprocess.run(
            [engine, "port", container, "8333/tcp"], env=engine_env, check=True, capture_output=True, text=True
        ).stdout.split()[0]
        endpoint = f"http://{port}"
        try:
            _create_buckets(endpoint)
        except (botocore.exceptions.BotoCoreError, botocore.exceptions.ClientError):
            logs = subprocess.run(
                [engine, "logs", container], env=engine_env, check=False, capture_output=True, text=True
            )
            print(f"$ {engine} logs {container}\n{logs.stdout}{logs.stderr}")
            raise
        yield endpoint
    finally:
        # the image declares /data a volume, which rm keeps without --volumes
        subprocess.run(
            [engine, "rm", "--force", "--volumes", container], env=engine_env, check=False, capture_output=True
        )


def _create_buckets(endpoint: str, timeout: float = 120) -> None:
    """Create the buckets once SeaweedFS takes writes, which is a few seconds after its start."""
    s3 = boto3.client("s3", endpoint_url=endpoint, aws_access_key_id=ADMIN[0], aws_secret_access_key=ADMIN[1])
    deadline = time.monotonic() + timeout
    while True:
        try:
            for bucket in BUCKETS:
                if bucket not in {b["Name"] for b in s3.list_buckets()["Buckets"]}:
                    s3.create_bucket(Bucket=bucket)
            s3.put_object(Bucket=INBOX, Key="ready", Body=b"")
            s3.delete_object(Bucket=INBOX, Key="ready")
            break
        except (botocore.exceptions.BotoCoreError, botocore.exceptions.ClientError):
            if time.monotonic() > deadline:
                raise
            time.sleep(1)
    for bucket in BUCKETS:
        s3.put_bucket_policy(Bucket=bucket, Policy=json.dumps(_bucket_policy(bucket)))


def _keys(keys: tuple[str, str]) -> dict[str, str]:
    return {"accessKey": keys[0], "secretKey": keys[1]}


def _bucket_policy(bucket: str) -> dict:
    """Grant each role its actions of :data:`PERMISSIONS` on *bucket*, and nothing else."""
    statements = []
    for role, buckets in PERMISSIONS.items():
        actions = buckets.get(bucket, ())
        principal = {"AWS": [f"arn:aws:iam::*:user/{credentials(role)[0]}"]}
        if object_actions := [a for a in actions if a != "s3:ListBucket"]:
            statements.append(
                {
                    "Effect": "Allow",
                    "Principal": principal,
                    "Action": object_actions,
                    "Resource": [f"arn:aws:s3:::{bucket}/*"],
                }
            )
        if "s3:ListBucket" in actions:
            statements.append(
                {
                    "Effect": "Allow",
                    "Principal": principal,
                    "Action": ["s3:ListBucket"],
                    "Resource": [f"arn:aws:s3:::{bucket}"],
                }
            )
    return {"Version": "2012-10-17", "Statement": statements}


def self_signed_certificate(directory: Path) -> tuple[Path, Path]:
    """Write a certificate for localhost and its key, for the BfArM mock.

    :returns: Paths of the certificate and of the key, both PEM.
    """
    key = ec.generate_private_key(ec.SECP256R1())
    name = x509.Name([x509.NameAttribute(NameOID.COMMON_NAME, "localhost")])
    now = datetime.datetime.now(datetime.UTC)
    certificate = (
        x509.CertificateBuilder()
        .subject_name(name)
        .issuer_name(name)
        .public_key(key.public_key())
        .serial_number(x509.random_serial_number())
        .not_valid_before(now - datetime.timedelta(minutes=5))
        .not_valid_after(now + datetime.timedelta(days=1))
        .add_extension(
            x509.SubjectAlternativeName([x509.DNSName("localhost"), x509.IPAddress(ipaddress.ip_address("127.0.0.1"))]),
            critical=False,
        )
        .add_extension(x509.BasicConstraints(ca=True, path_length=None), critical=True)
        .sign(key, hashes.SHA256())
    )
    cert_path, key_path = directory / "bfarm-mock.crt", directory / "bfarm-mock.key"
    cert_path.write_bytes(certificate.public_bytes(serialization.Encoding.PEM))
    key_path.write_bytes(
        key.private_bytes(serialization.Encoding.PEM, serialization.PrivateFormat.PKCS8, serialization.NoEncryption())
    )
    return cert_path, key_path


class BfarmMock:
    """HTTPS stand-in for the BfArM Prüfbericht API.

    It hands out a token for any client credentials. It answers each Prüfbericht with 200 while ``accept``
    is true, else with 400, and keeps every accepted one in ``received``.
    """

    def __init__(self, cert_path: Path, key_path: Path):
        self.accept = True
        self.received: list[dict] = []
        mock = self

        class Handler(BaseHTTPRequestHandler):
            def do_POST(self):
                body = self.rfile.read(int(self.headers.get("Content-Length", 0)))
                if self.path == "/token":
                    self._reply(200, {"access_token": "e2e-token", "expires_in": 3600})
                elif self.path == "/api/upload" and mock.accept:
                    mock.received.append(json.loads(body))
                    self._reply(200, {"status": "success"})
                elif self.path == "/api/upload":
                    self._reply(400, {"status": "rejected"})
                else:
                    self._reply(404, {})

            def _reply(self, status: int, payload: dict) -> None:
                data = json.dumps(payload).encode()
                self.send_response(status)
                self.send_header("Content-Type", "application/json")
                self.send_header("Content-Length", str(len(data)))
                self.end_headers()
                self.wfile.write(data)

            def log_message(self, format, *args):
                pass

        context = ssl.SSLContext(ssl.PROTOCOL_TLS_SERVER)
        context.load_cert_chain(cert_path, key_path)
        self._server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
        self._server.socket = context.wrap_socket(self._server.socket, server_side=True)
        self.url = f"https://127.0.0.1:{self._server.server_address[1]}"
        self._thread = threading.Thread(target=self._server.serve_forever, daemon=True)
        self._thread.start()

    def reset(self) -> None:
        self.accept = True
        self.received.clear()

    def close(self) -> None:
        self._server.shutdown()
        self._server.server_close()
