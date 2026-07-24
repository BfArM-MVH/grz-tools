"""
End-to-end integration test for `grzctl process` with validation ENABLED.

The existing happy-path integration tests in ``test_process.py`` all run with
``--no-validate``. This module fills the gap: it drives a *valid* pretend
submission through the complete streaming pipeline

    inbox download -> decrypt -> validate (fastq/bam/checksum) -> re-encrypt -> archive

and then proves correctness by decrypting the re-encrypted archive output and
comparing it byte-for-byte (via SHA256) against the original plaintext files.

This is an "on the fly" test added during review; it does not modify any
branch code.
"""

from __future__ import annotations

import hashlib
import io
from pathlib import Path

import boto3
import click.testing
import crypt4gh.keys
import crypt4gh.lib
import grzctl.cli
import pytest
import yaml
from moto import mock_aws

MOCK_FILES_DIR = Path(__file__).parent.parent / "mock_files"
VALID_SUBMISSION_DIR = MOCK_FILES_DIR / "submissions" / "valid_submission"
SUBMISSION_ID = "260914050_2024-07-15_c64603a7"


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _config_content(
    grz_private_key_path: Path,
    grz_public_key_path: Path,
    db_alice_private_key_file_path: Path,
    db_known_keys_file_path: Path,
    db_file: Path,
) -> dict:
    """Minimal config: one inbox + consented/non-consented archives sharing the GRZ keypair."""
    s3_common = {
        "endpoint_url": "https://s3.amazonaws.com",
        "access_key": "testing",
        "secret": "testing",
    }
    return {
        "s3": {**s3_common, "bucket": "inbox"},
        "consented_archive_s3": {**s3_common, "bucket": "consented-archive"},
        "non_consented_archive_s3": {**s3_common, "bucket": "non-consented-archive"},
        "keys": {
            "grz_private_key_path": str(grz_private_key_path),
            "consented_archive_public_key_path": str(grz_public_key_path),
            "non_consented_archive_public_key_path": str(grz_public_key_path),
        },
        "pruefbericht": {
            "authorization_url": "https://bfarm.localhost/token",
            "api_base_url": "https://bfarm.localhost/api/",
            "client_id": "pytest",
            "client_secret": "pysecret",
        },
        "db": {
            "database_url": f"sqlite:///{db_file}",
            "author": {"name": "Alice", "private_key_path": str(db_alice_private_key_file_path)},
            "known_public_keys": str(db_known_keys_file_path),
        },
    }


def _upload_submission_to_inbox(inbox_bucket) -> None:
    inbox_bucket.upload_file(
        Filename=str(VALID_SUBMISSION_DIR / "metadata" / "metadata.json"),
        Key=f"{SUBMISSION_ID}/metadata/metadata.json",
    )
    for encrypted_file in (VALID_SUBMISSION_DIR / "encrypted_files").glob("*.c4gh"):
        inbox_bucket.upload_file(
            Filename=str(encrypted_file),
            Key=f"{SUBMISSION_ID}/files/{encrypted_file.name}",
        )


def _decrypt_c4gh_bytes(encrypted: bytes, private_key: bytes) -> bytes:
    """Decrypt crypt4gh bytes using the reference crypt4gh library (not the streaming code)."""
    out = io.BytesIO()
    crypt4gh.lib.decrypt([(0, private_key, None)], io.BytesIO(encrypted), out, sender_pubkey=None)
    return out.getvalue()


def _encrypt_c4gh_bytes(plaintext: bytes, recipient_public_key: bytes) -> bytes:
    """Encrypt bytes for a recipient using the reference crypt4gh library."""
    from nacl.public import PrivateKey

    sender_sk = bytes(PrivateKey.generate())
    out = io.BytesIO()
    crypt4gh.lib.encrypt([(0, sender_sk, recipient_public_key)], io.BytesIO(plaintext), out)
    return out.getvalue()


@pytest.fixture
def _aws(monkeypatch):
    monkeypatch.setenv("AWS_ACCESS_KEY_ID", "testing")
    monkeypatch.setenv("AWS_SECRET_ACCESS_KEY", "testing")
    monkeypatch.setenv("MOTO_ALLOW_NONEXISTENT_REGION", "1")
    with mock_aws():
        yield


@pytest.fixture
def buckets(_aws):
    conn = boto3.client("s3")
    for name in ("inbox", "consented-archive", "non-consented-archive"):
        conn.create_bucket(Bucket=name)
    s3 = boto3.resource("s3")
    return {
        "inbox": s3.Bucket("inbox"),
        "consented": s3.Bucket("consented-archive"),
        "non_consented": s3.Bucket("non-consented-archive"),
    }


@pytest.fixture
def config_file(
    tmp_path,
    crypt4gh_grz_private_key_file_path,
    crypt4gh_grz_public_key_file_path,
    db_alice_private_key_file_path,
    db_known_keys_file_path,
):
    db_file = tmp_path / "db" / "test.db"
    db_file.parent.mkdir(parents=True, exist_ok=True)
    content = _config_content(
        crypt4gh_grz_private_key_file_path,
        crypt4gh_grz_public_key_file_path,
        db_alice_private_key_file_path,
        db_known_keys_file_path,
        db_file,
    )
    path = tmp_path / "config.process.yaml"
    path.write_text(yaml.dump(content))
    return path


def _run_process(config_file: Path, output_dir: Path, *extra: str, catch_exceptions: bool = False):
    runner = click.testing.CliRunner()
    cli = grzctl.cli.build_cli()
    args = [
        "process",
        "--config-file",
        str(config_file),
        "--submission-id",
        SUBMISSION_ID,
        "--output-dir",
        str(output_dir),
        "--no-submit-pruefbericht",
        "--no-update-db",
        *extra,
    ]
    return runner.invoke(cli, args, catch_exceptions=catch_exceptions)


class TestProcessWithValidationEndToEnd:
    def test_valid_submission_passes_validation_and_archives_correct_content(
        self,
        buckets,
        config_file,
        crypt4gh_grz_private_key_file_path,
        tmp_path,
    ):
        """
        Drive a valid submission through the full pipeline WITH validation enabled,
        then verify the archived (re-encrypted) files decrypt back to the originals.
        """
        _upload_submission_to_inbox(buckets["inbox"])

        result = _run_process(config_file, tmp_path / "work", "--validate")
        assert result.exit_code == 0, f"process --validate failed:\n{result.output}"

        # --- archive landed in consented bucket (valid_submission has research consent) ---
        consented_keys = {o.key for o in buckets["consented"].objects.all()}
        assert any(k.endswith("metadata/metadata.json") for k in consented_keys), consented_keys
        archived_c4gh = [k for k in consented_keys if k.endswith(".c4gh")]
        assert archived_c4gh, f"no encrypted files archived: {consented_keys}"

        # nothing leaked into the non-consented archive
        assert {o.key for o in buckets["non_consented"].objects.all()} == set()

        # --- prove round-trip correctness: archived ciphertext -> plaintext == original ---
        grz_priv = crypt4gh.keys.get_private_key(str(crypt4gh_grz_private_key_file_path), lambda: "")
        s3 = boto3.client("s3")

        original_checksums = {
            p.name: _sha256(p.read_bytes())
            for p in (VALID_SUBMISSION_DIR / "files").glob("*")
            if p.is_file()
        }

        verified = 0
        for key in archived_c4gh:
            body = s3.get_object(Bucket="consented-archive", Key=key)["Body"].read()
            decrypted = _decrypt_c4gh_bytes(body, grz_priv)
            original_name = Path(key).name[: -len(".c4gh")]
            assert original_name in original_checksums, f"unexpected archived file {original_name}"
            assert _sha256(decrypted) == original_checksums[original_name], (
                f"re-encrypted archive content for {original_name} does not match original after decrypt"
            )
            verified += 1

        assert verified == len(archived_c4gh) > 0
        # the submission contains gzipped FASTQ -> the decompress+FastqValidator path was exercised
        assert any(name.endswith(".fastq.gz") for name in original_checksums)

    def test_validation_records_per_file_success_state(
        self,
        buckets,
        config_file,
        tmp_path,
    ):
        """The processing progress log should mark every file as successfully processed."""
        import json as _json

        _upload_submission_to_inbox(buckets["inbox"])
        work = tmp_path / "work"
        result = _run_process(config_file, work, "--validate")
        assert result.exit_code == 0, result.output

        progress = work / "logs" / "progress_processing.cjson"
        assert progress.is_file(), "expected processing progress log"

        states = []
        for line in progress.read_text().splitlines():
            line = line.strip()
            if line:
                states.append(_json.loads(line))
        assert states, "progress log is empty"

        def _is_success(entry):
            state = entry.get("state", entry)
            return bool(state.get("processing_successful"))

        assert all(_is_success(s) for s in states), f"some files not marked successful: {states}"
        # all submission files (fastq/vcf/bed) were processed
        assert len(states) >= 6


class TestProcessValidationRejection:
    """A corrupted (but properly encrypted) file must fail validation and abort the upload."""

    def test_corrupted_file_fails_validation_and_archives_nothing(
        self,
        buckets,
        config_file,
        crypt4gh_grz_public_key_file_path,
        tmp_path,
    ):
        _upload_submission_to_inbox(buckets["inbox"])

        # pick one real, encrypted FASTQ inbox object and overwrite it with
        # properly-encrypted-but-corrupted content (3 lines -> not a multiple of 4,
        # and a checksum/size that won't match the metadata).
        s3 = boto3.client("s3")
        target_key = next(
            o.key for o in buckets["inbox"].objects.all() if o.key.endswith(".fastq.gz.c4gh")
        )

        import gzip

        bad_plaintext = gzip.compress(b"@read1\nACGT\n+\n")  # 3 lines, wrong content
        recipient_pub = crypt4gh.keys.get_public_key(str(crypt4gh_grz_public_key_file_path))
        s3.put_object(
            Bucket="inbox",
            Key=target_key,
            Body=_encrypt_c4gh_bytes(bad_plaintext, recipient_pub),
        )

        # use the real-CLI default (exceptions caught -> non-zero exit) rather than
        # re-raising, so we exercise the command's actual failure exit behavior.
        result = _run_process(config_file, tmp_path / "work", "--validate", catch_exceptions=True)

        # pipeline must fail (validation rejected the corrupted FASTQ)
        assert result.exit_code != 0, f"expected validation failure, got success:\n{result.output}"
        assert isinstance(result.exception, BaseException)

        # and nothing may remain in either archive (incomplete submission cleaned up)
        assert {o.key for o in buckets["consented"].objects.all()} == set()
        assert {o.key for o in buckets["non_consented"].objects.all()} == set()
