"""The LE's and the GRZ's side of one end-to-end scenario."""

import json
import os
import secrets
import shutil
import subprocess
import sys
from collections import defaultdict
from collections.abc import Callable
from pathlib import Path

import boto3
import yaml
from grz_common.workers.submission import SubmissionMetadata
from grz_db.models.submission import SubmissionDb

from . import services
from .scenarios import STEPS, Flow, Run, Step, Submission

MOCK_FILES = Path(__file__).parents[1] / "mock_files"
GRZ_ID = "GRZM00123"
LE_ID = "123456789"
"""The IDs in the metadata of grz-mini-test-data."""

PART_SIZE = 5 * 1024**2
"""The smallest part size S3 allows, so that the files of a grown mini submission go up in several parts."""

_BIN = Path(sys.executable).parent
_QC_OK = (
    "mkdir -p '{output_basepath}' && find '{submission_basepath}/files' -type f | sort > '{output_basepath}/files.txt'"
)
_QC_FAILS = "exit 3"
TRANSITION_WARNINGS = ("Expected any of", "does not contain any of")
"""Parts of the warnings that grzctl logs for a state transition it did not expect."""


class Grz:
    """One GRZ with its own database and QC storage, on the shared buckets, and its LE.

    Each submission gets a random tanG, so scenarios never meet in the buckets.
    """

    def __init__(
        self,
        work_dir: Path,
        s3_endpoint: str,
        database_url: str,
        bfarm: services.BfarmMock,
        ca_bundle: Path,
        datasets: Callable[[str], Path],
    ):
        self.work_dir = work_dir
        self._datasets = datasets
        self.s3_endpoint = s3_endpoint
        self.database_url = database_url
        self.bfarm = bfarm
        self.qc_storage = work_dir / "qc"
        self.outputs: list[str] = []
        self._tans: defaultdict[str, str] = defaultdict(lambda: secrets.token_hex(32))
        self._env = {k: v for k, v in os.environ.items() if not k.lower().startswith(("grz_", "c4gh_", "aws_"))}
        self._env["REQUESTS_CA_BUNDLE"] = str(ca_bundle)

    # LE

    def upload(self, submission: Submission, name: str) -> str:
        """Prepare *submission* from its mini dataset and upload it with grz-cli, as the LE does.

        :param name: Directory name of this submission in the work directory.
        :returns: The submission ID.
        """
        source = self._datasets(submission.dataset)
        submission_dir = self.work_dir / name
        shutil.copytree(source, submission_dir)
        metadata_path = submission_dir / "metadata" / "metadata.json"
        metadata = json.loads(metadata_path.read_text())
        metadata["submission"]["tanG"] = self._tans[submission.tan] if submission.tan else secrets.token_hex(32)
        for change in submission.metadata:
            change(metadata)
        for file_change in submission.files:
            file_change(submission_dir / "files", metadata)
        metadata_path.write_text(json.dumps(metadata, indent=2))

        config = self._write_config(
            f"grz-cli-{name}.yaml",
            {
                "identifiers": {"grz": GRZ_ID, "le": LE_ID},
                "keys": {
                    "grz_public_key_path": str(
                        MOCK_FILES
                        / ("grz_mock_public_key.pub" if submission.encrypt_for_grz else "submitter_mock_public_key.pub")
                    )
                },
                "s3": {
                    "endpoint_url": self.s3_endpoint,
                    "bucket": services.INBOX,
                    "access_key": services.credentials("le")[0],
                    "secret": services.credentials("le")[1],
                    "multipart_chunksize": PART_SIZE,
                },
            },
        )
        common = ["--config-file", str(config), "--submission-dir", str(submission_dir)]
        if submission.client_validation:
            self.run("grz-cli", "submit", *common, "--threads", "2")
        else:
            self.run("grz-cli", "encrypt", *common, "--no-check-validation-logs")
            self.run("grz-cli", "upload", *common, "--threads", "2")

        submission_id = SubmissionMetadata(metadata_path).content.submission_id
        for path in submission.removed_from_inbox:
            self.s3().Object(services.INBOX, f"{submission_id}/files/{path}.c4gh").delete()
        return submission_id

    # GRZ

    def run_flow(self, flow: Flow, submission_id: str, runs: tuple[Run, ...]) -> None:
        """Run the GRZ's *flow* over the submission once per run; a rerun goes on where the previous run failed."""
        done: set[str] = set()
        for index, run in enumerate(runs):
            self.bfarm.accept = run.bfarm_accepts
            if flow == "process":
                self._process(submission_id, run)
            else:
                self._steps(submission_id, run, done, first=index == 0)

    def _process(self, submission_id: str, run: Run) -> None:
        """Run ``grzctl process``; a rerun reuses the output directory and its progress logs."""
        self.run(
            "grzctl",
            "--config",
            str(self._grzctl_config(run)),
            "process",
            "--submission-id",
            submission_id,
            "--output-dir",
            str(self.work_dir / "process" / submission_id),
            "--threads",
            "2",
            check=run.fails_at is None,
        )

    def _steps(self, submission_id: str, run: Run, done: set[str], first: bool) -> None:
        """Run the single commands of the steps not yet *done*, up to the one that fails."""
        config = ["--config", str(self._grzctl_config(run))]
        submission_dir = self.work_dir / "steps" / submission_id
        in_dir = ["--submission-dir", str(submission_dir), "--update-db"]
        if first:
            submission_dir.parent.mkdir(exist_ok=True)
            self.run("grzctl", *config, "db", "sync-from-inbox", "--submitter-id", LE_ID, "--inbox", services.INBOX)
        commands: dict[Step, Callable[[bool], None]] = {
            "download": lambda ok: self.run(
                "grzctl", *config, "download", "--submission-id", submission_id, "--output-dir", str(submission_dir),
                "--update-db", "--populate", check=ok,
            ),
            "decrypt": lambda ok: self.run("grzctl", *config, "decrypt", *in_dir, check=ok),
            "validate": lambda ok: self.run("grzctl", *config, "validate", *in_dir, "--submitter-id", LE_ID, check=ok),
            # the re-encrypted files replace the downloaded ones in encrypted_files/
            "encrypt": lambda ok: self.run("grzctl", *config, "encrypt", *in_dir, "--force", check=ok),
            "archive": lambda ok: self.run("grzctl", *config, "archive", *in_dir, check=ok),
            "report": lambda ok: self._report(config, submission_id, ok),
            "detailed_qc": lambda ok: self._detailed_qc(config, submission_id, run, ok),
            "clean": lambda ok: self.run(
                "grzctl", *config, "clean", "--submission-id", submission_id, "--yes-i-really-mean-it", "--update-db", check=ok
            ),
        }  # fmt: skip
        for step in STEPS:
            if step in done:
                continue
            commands[step](run.fails_at != step)
            if run.fails_at == step:
                return
            done.add(step)

    def _report(self, config: list[str], submission_id: str, ok: bool) -> None:
        pruefbericht = self.run(
            "grzctl", *config, "pruefbericht", "generate", "from-database", "--submission-id", submission_id
        )
        path = self.work_dir / f"pruefbericht-{submission_id}.json"
        path.write_text(pruefbericht.stdout)
        self.run(
            "grzctl", *config, "pruefbericht", "submit", "--pruefbericht-file", str(path), "--submission-id", submission_id,
            "--update-db", check=ok,
        )  # fmt: skip

    def _detailed_qc(self, config: list[str], submission_id: str, run: Run, ok: bool) -> None:
        """Decide and run detailed QC by hand, as a GDC does without ``grzctl process``."""
        decision = self.run("grzctl", *config, "db", "should-qc", submission_id)
        if decision.stdout.strip() != "true":
            return
        update = ["grzctl", *config, "db", "submission", "update", submission_id]
        self.run(*update, "qcing", "--ignore-error-state")
        submission_dir = self.work_dir / "steps" / submission_id
        command = (_QC_OK if run.qc_succeeds else _QC_FAILS).format(
            submission_basepath=submission_dir, output_basepath=self.qc_output("steps", submission_id).parent
        )
        succeeded = subprocess.run(command, shell=True, check=False).returncode == 0
        if succeeded:
            self.run(*update, "qced")
        else:
            self.run(*update, "error", "--failure-reason", "detailed_qc_error")
        assert succeeded == ok, f"QC workflow {'succeeded' if succeeded else 'failed'}"

    def qc_output(self, flow: Flow, submission_id: str) -> Path:
        """Return the file in which the QC workflow lists the files it got."""
        if flow == "process":
            return self.qc_storage / submission_id / "qc" / "files.txt"
        return self.work_dir / "steps" / submission_id / "qc" / "files.txt"

    def run(self, command: str, *args: str, check: bool = True) -> subprocess.CompletedProcess:
        """Run grz-cli or grzctl of this environment and keep its output for the test report."""
        result = subprocess.run(
            [str(_BIN / command), *args], env=self._env, capture_output=True, text=True, check=False, timeout=600
        )
        self.outputs.append(f"$ {command} {' '.join(args)}\n[exit {result.returncode}]\n{result.stdout}{result.stderr}")
        print(self.outputs[-1])
        if (result.returncode == 0) != check:
            raise AssertionError(
                f"{command} {args[:4]} exited with {result.returncode}, expected {'0' if check else 'non-zero'}"
            )
        return result

    # inspection

    def states(self, submission_id: str) -> tuple[list[str], str | None]:
        """Return the submission's state history and the failure reason of its last state."""
        submission = SubmissionDb(db_url=self.database_url, author=None).get_submission(submission_id)
        assert submission is not None, f"database holds no submission {submission_id}"
        logs = sorted(submission.states, key=lambda log: (log.timestamp, log.id))
        reason = logs[-1].failure_reason
        return [log.state.value for log in logs], reason.value if reason else None

    def submission(self, submission_id: str):
        return SubmissionDb(db_url=self.database_url, author=None).get_submission(submission_id)

    def s3(self):
        return boto3.resource(
            "s3",
            endpoint_url=self.s3_endpoint,
            aws_access_key_id=services.ADMIN[0],
            aws_secret_access_key=services.ADMIN[1],
        )

    def keys(self, bucket: str, submission_id: str) -> set[str]:
        return {obj.key for obj in self.s3().Bucket(bucket).objects.filter(Prefix=f"{submission_id}/")}

    # config

    def _grzctl_config(self, run: Run) -> Path:
        archive = lambda name, bucket: {
            "s3": {**self._s3(name), "bucket": bucket},
            "public_key_path": str(MOCK_FILES / f"archive_{name}.pub"),
            "private_key_path": str(MOCK_FILES / f"archive_{name}.sec"),
        }
        return self._write_config(
            "grzctl.yaml",
            {
                "leistungserbringer": {
                    LE_ID: {
                        "inbox_buckets": {
                            services.INBOX: {
                                **self._s3("inbox"),
                                "private_key_path": str(MOCK_FILES / "grz_mock_private_key.sec"),
                            }
                        }
                    }
                },
                "archives": {
                    "consented": archive("consented", services.ARCHIVE_CONSENTED),
                    "non_consented": archive("non_consented", services.ARCHIVE_NON_CONSENTED),
                    "interrogation": {"s3": {**self._s3("interrogation"), "bucket": services.ARCHIVE_INTERROGATION}},
                },
                "pruefbericht": {
                    "authorization_url": f"{self.bfarm.url}/token",
                    "client_id": "e2e",
                    "client_secret": "e2e",
                    "api_base_url": f"{self.bfarm.url}/api",
                },
                "detailed_qc": {
                    "local_storage": str(self.qc_storage),
                    "salt": "e2e",
                    "target_percentage": 100,
                    "auto_run": True,
                    "shell_command": _QC_OK if run.qc_succeeds else _QC_FAILS,
                },
                "db": {
                    "database_url": self.database_url,
                    "author": {
                        "name": "Alice",
                        "private_key_path": str(MOCK_FILES / "db" / "alice_mock_private_key.sec"),
                        "private_key_passphrase": "test",
                    },
                    "known_public_keys_file": str(MOCK_FILES / "db" / "known_keys"),
                },
                "identifiers": {"grz": GRZ_ID},
            },
        )

    def _s3(self, role: str) -> dict[str, str]:
        """The S3 connection of *role*, whose bucket policies grant exactly the actions of ``PERMISSIONS``."""
        access_key, secret = services.credentials(role)
        return {
            "endpoint_url": self.s3_endpoint,
            "access_key": access_key,
            "secret": secret,
            "multipart_chunksize": PART_SIZE,
        }

    def _write_config(self, name: str, content: dict) -> Path:
        path = self.work_dir / name
        path.write_text(yaml.safe_dump(content))
        return path
