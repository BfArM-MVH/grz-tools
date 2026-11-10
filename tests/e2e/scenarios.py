"""End-to-end scenarios: one submission from grz-mini-test-data, how the LE uploads it, and what the GRZ ends with.

Each scenario changes a copy of a mini submission, uploads it with grz-cli, and lets the GRZ process it in each of
its flows: ``grzctl process``, or the single commands from ``download`` to ``clean`` that GDCs run by hand or in
their own automation. A scenario says where a run fails, not which states it records; :func:`expected_states`
derives those per flow. A scenario's description says what it shows; its fields say how.
"""

import gzip
import hashlib
import importlib.resources
import json
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

from grz_pydantic_models_testing import example_research_consent

MetadataChange = Callable[[dict], None]
FileChange = Callable[[Path, dict], None]
"""Changes a file in ``files/`` of the submission; gets the directory and the metadata to keep it consistent."""

Flow = Literal["process", "steps"]
Step = Literal["download", "decrypt", "validate", "encrypt", "archive", "report", "detailed_qc", "clean"]
STEPS: tuple[Step, ...] = ("download", "decrypt", "validate", "encrypt", "archive", "report", "detailed_qc", "clean")
"""The steps of the GRZ in order. ``grzctl process`` does the first five in one streaming pass."""
_ARCHIVAL: tuple[Step, ...] = STEPS[:5]
_STATES: dict[str, tuple[str, str]] = {
    "process": ("Processing", "Processed"),
    "download": ("Downloading", "Downloaded"),
    "decrypt": ("Decrypting", "Decrypted"),
    "validate": ("Validating", "Validated"),
    "encrypt": ("Encrypting", "Encrypted"),
    "archive": ("Archiving", "Archived"),
    "report": ("Reporting", "Reported"),
    "detailed_qc": ("QCing", "QCed"),
    "clean": ("Cleaning", "Cleaned"),
}


def _files(metadata: dict) -> list[dict]:
    return [f for donor in metadata["donors"] for lab in donor["labData"] for f in lab["sequenceData"]["files"]]


def _file_entry(metadata: dict, path: str) -> dict:
    return next(f for f in _files(metadata) if f["filePath"] == path)


def initial(metadata: dict) -> None:
    """Make it an initial submission."""
    metadata["submission"]["submissionType"] = "initial"


def research_consent(metadata: dict) -> None:
    """Give every donor a research consent that permits research at the submission date."""
    scope = json.loads(
        importlib.resources.files(example_research_consent).joinpath("minimal_consented.json").read_text()
    )
    for donor in metadata["donors"]:
        donor["researchConsents"][0]["scope"] = scope


def submission_date(date: str) -> MetadataChange:
    """Move the submission to another date, which gives it another submission ID."""

    def change(metadata: dict) -> None:
        metadata["submission"]["submissionDate"] = date

    change.__name__ = f"submission_date({date})"
    return change


def wrong_checksum(path: str) -> MetadataChange:
    """State a checksum in metadata.json that does not match the file."""

    def change(metadata: dict) -> None:
        _file_entry(metadata, path)["fileChecksum"] = hashlib.sha256(b"not the file").hexdigest()

    change.__name__ = f"wrong_checksum({path})"
    return change


def _restate(files_dir: Path, metadata: dict, path: str) -> None:
    """Update checksum and size in metadata.json to the changed file, so only the format check can catch it."""
    data = (files_dir / path).read_bytes()
    entry = _file_entry(metadata, path)
    entry["fileChecksum"] = hashlib.sha256(data).hexdigest()
    entry["fileSizeInBytes"] = len(data)


def cut_mid_record(path: str) -> FileChange:
    """Cut a gzipped FASTQ in the middle of a record, as a broken copy would."""

    def change(files_dir: Path, metadata: dict) -> None:
        lines = gzip.decompress((files_dir / path).read_bytes()).splitlines(keepends=True)
        (files_dir / path).write_bytes(gzip.compress(b"".join(lines[: len(lines) // 2 + 2])))
        _restate(files_dir, metadata, path)

    change.__name__ = f"cut_mid_record({path})"
    return change


def drop_reads(path: str, count: int) -> FileChange:
    """Drop the last *count* reads of a gzipped FASTQ, so its mate has more reads."""

    def change(files_dir: Path, metadata: dict) -> None:
        lines = gzip.decompress((files_dir / path).read_bytes()).splitlines(keepends=True)
        (files_dir / path).write_bytes(gzip.compress(b"".join(lines[: -4 * count])))
        _restate(files_dir, metadata, path)

    change.__name__ = f"drop_reads({path}, {count})"
    return change


def grow_beyond_parts(*paths: str, copies: int = 3) -> FileChange:
    """Repeat the reads of gzipped FASTQs *copies* times, so each file spans several parts of a multipart upload.

    Real submission files are many GB, but the mini ones fit into one part. Each copy gets its own read names,
    and mates keep matching names.
    """

    def change(files_dir: Path, metadata: dict) -> None:
        for path in paths:
            lines = gzip.decompress((files_dir / path).read_bytes()).splitlines(keepends=True)
            grown = []
            for copy in range(copies):
                for index, line in enumerate(lines):
                    if index % 4 == 0:
                        name, _, mate = line.rstrip(b"\n").rpartition(b"/")
                        line = b"%s_%d/%s\n" % (name, copy, mate)
                    grown.append(line)
            (files_dir / path).write_bytes(gzip.compress(b"".join(grown), compresslevel=1, mtime=0))
            _restate(files_dir, metadata, path)

    change.__name__ = f"grow_beyond_parts({', '.join(paths)})"
    return change


@dataclass(frozen=True)
class Submission:
    """A mini submission, changed and uploaded by the LE."""

    dataset: Literal["panel", "wgs", "wgs_lr"]
    metadata: tuple[MetadataChange, ...] = ()
    files: tuple[FileChange, ...] = ()
    tan: str | None = None
    """Submissions of a scenario with the same label share one random tanG; without one, a submission gets its own."""
    client_validation: bool = True
    """False: the LE skips ``grz-cli validate`` and encrypts with ``--no-check-validation-logs``."""
    encrypt_for_grz: bool = True
    """False: the LE encrypts for another key than the GRZ's."""
    removed_from_inbox: tuple[str, ...] = ()
    """Files (``filePath``) deleted from the inbox after the upload."""


@dataclass(frozen=True)
class Run:
    """One run of the GRZ over the submission. A rerun goes on where the previous run failed."""

    fails_at: Step | None = None
    """The step at which this run stops with an error; ``None`` if it gets through."""
    bfarm_accepts: bool = True
    qc_succeeds: bool = True


@dataclass(frozen=True)
class Outcome:
    """What the GRZ holds after the last run."""

    failure_reason: str | None = None
    """Failure reason of the last state, for a submission that ends in ``Error``."""
    archive: Literal["consented", "non_consented"] | None = None
    """Archive that holds the submission; ``None`` if neither does."""
    inbox_cleaned: bool = False
    pruefberichte: int = 0
    """Prüfberichte that BfArM accepted."""
    detailed_qc: bool = False
    """Whether the submission is selected for detailed QC (all initial ones, at the target of 100 %)."""


@dataclass(frozen=True)
class Scenario:
    id: str
    description: str
    submission: Submission
    outcome: Outcome
    runs: tuple[Run, ...] = (Run(),)
    earlier: tuple[Submission, ...] = field(default=())
    """Submissions that the GRZ processed successfully before, in the same flow and database."""


def expected_states(scenario: Scenario, flow: Flow) -> list[str]:
    """Return the state history that *flow* records for *scenario*.

    ``process`` records one state pair for the archival steps, and a rerun records it again without streaming.
    The single commands start from ``Uploaded``, which ``grzctl db sync-from-inbox`` records, and a rerun starts
    at the step that failed.
    """
    steps = [step for step in STEPS if step != "detailed_qc" or scenario.outcome.detailed_qc]
    states = ["Uploaded"] if flow == "steps" else []
    done: set[str] = set()
    for run in scenario.runs:
        if flow == "process":
            pending = ["process", *(step for step in steps[len(_ARCHIVAL) :] if step not in done)]
        else:
            pending = [step for step in steps if step not in done]
        for step in pending:
            start, end = _STATES[step]
            if run.fails_at == step or (step == "process" and run.fails_at in _ARCHIVAL):
                states += [start, "Error"]
                break
            states += [start, end]
            done.update(_ARCHIVAL if step == "process" else {step})
    return states


SCENARIOS = (
    Scenario(
        id="panel-test",
        description="A test submission without research consent is archived non-consented, reported and cleaned. "
        "Only initial submissions take part in detailed QC.",
        submission=Submission("panel"),
        outcome=Outcome(archive="non_consented", inbox_cleaned=True, pruefberichte=1),
    ),
    Scenario(
        id="wgs-trio-initial",
        description="An initial trio submission with research consent of all donors is archived consented. "
        "The QC target of 100 % selects it, and the QC workflow runs after the Prüfbericht.",
        submission=Submission("wgs", metadata=(initial, research_consent)),
        outcome=Outcome(archive="consented", inbox_cleaned=True, pruefberichte=1, detailed_qc=True),
    ),
    Scenario(
        id="wgs-lr-test",
        description="A long-read submission with BAM files goes through like a short-read one.",
        submission=Submission("wgs_lr"),
        outcome=Outcome(archive="non_consented", inbox_cleaned=True, pruefberichte=1),
    ),
    Scenario(
        id="checksum-mismatch",
        description="metadata.json states a wrong checksum, and the LE skipped validation. "
        "Nothing is archived or reported, and the inbox keeps the submission.",
        submission=Submission("panel", metadata=(wrong_checksum("tumor_1.fq.gz"),), client_validation=False),
        runs=(Run(fails_at="validate"),),
        outcome=Outcome(failure_reason="validation_error"),
    ),
    Scenario(
        id="fastq-cut-mid-record",
        description="A FASTQ ends in the middle of a record. Checksum and size match the broken file, "
        "so only the FASTQ check catches it.",
        submission=Submission("panel", files=(cut_mid_record("tumor_1.fq.gz"),), client_validation=False),
        runs=(Run(fails_at="validate"),),
        outcome=Outcome(failure_reason="validation_error"),
    ),
    Scenario(
        id="mates-differ-in-read-count",
        description="R2 of a pair has fewer reads than R1.",
        submission=Submission("panel", files=(drop_reads("tumor_2.fq.gz", 10),), client_validation=False),
        runs=(Run(fails_at="validate"),),
        outcome=Outcome(failure_reason="validation_error"),
    ),
    Scenario(
        id="encrypted-for-another-key",
        description="The LE encrypted the files for another key than the GRZ's.",
        submission=Submission("panel", encrypt_for_grz=False),
        runs=(Run(fails_at="decrypt"),),
        outcome=Outcome(failure_reason="decryption_error"),
    ),
    Scenario(
        id="file-missing-in-inbox",
        description="A file that metadata.json lists is missing in the inbox.",
        submission=Submission("panel", removed_from_inbox=("tumor_2.fq.gz",)),
        runs=(Run(fails_at="download"),),
        outcome=Outcome(failure_reason="file_not_found"),
    ),
    Scenario(
        id="duplicate-tang",
        description="An earlier submission already used the tanG.",
        earlier=(Submission("panel", tan="shared"),),
        submission=Submission("panel", tan="shared", metadata=(submission_date("2024-11-09"),)),
        runs=(Run(fails_at="download"),),
        outcome=Outcome(failure_reason="duplicate_tang"),
    ),
    Scenario(
        id="duplicate-initial",
        description="The case already has an initial submission that passed basic QC.",
        earlier=(Submission("panel", metadata=(initial,)),),
        submission=Submission("panel", metadata=(initial,)),
        runs=(Run(fails_at="validate"),),
        outcome=Outcome(failure_reason="duplicate_initial"),
    ),
    Scenario(
        id="pruefbericht-rejected-then-rerun",
        description="BfArM rejects the Prüfbericht, so the inbox keeps the archived submission. "
        "A rerun after BfArM accepts reports and cleans it without archiving it again.",
        submission=Submission("panel"),
        runs=(Run(fails_at="report", bfarm_accepts=False), Run()),
        outcome=Outcome(archive="non_consented", inbox_cleaned=True, pruefberichte=1),
    ),
    Scenario(
        id="files-in-parts",
        description="FASTQs larger than one part go through multipart uploads, and a multipart copy into the archive.",
        submission=Submission("panel", files=(grow_beyond_parts("tumor_1.fq.gz", "tumor_2.fq.gz"),)),
        outcome=Outcome(archive="non_consented", inbox_cleaned=True, pruefberichte=1),
    ),
    Scenario(
        id="file-in-parts-cut-mid-record",
        description="A FASTQ larger than one part ends in the middle of a record. "
        "grzctl process aborts its multipart upload to the interrogation bucket, so no parts stay behind.",
        submission=Submission(
            "panel",
            files=(grow_beyond_parts("tumor_1.fq.gz", "tumor_2.fq.gz"), cut_mid_record("tumor_1.fq.gz")),
            client_validation=False,
        ),
        runs=(Run(fails_at="validate"),),
        outcome=Outcome(failure_reason="validation_error"),
    ),
    Scenario(
        id="detailed-qc-fails-then-rerun",
        description="The QC workflow fails, so the inbox keeps the submission. A rerun runs it again and cleans.",
        submission=Submission("panel", metadata=(initial,)),
        runs=(Run(fails_at="detailed_qc", qc_succeeds=False), Run()),
        outcome=Outcome(archive="non_consented", inbox_cleaned=True, pruefberichte=1, detailed_qc=True),
    ),
)
