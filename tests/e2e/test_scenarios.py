"""Run every scenario of :mod:`scenarios` in both flows of the GRZ, from the LE's upload to the GRZ's last step."""

import hashlib
import io

import crypt4gh.keys
import crypt4gh.lib
import pytest
from grz_common.workers.submission import SubmissionMetadata

from . import services
from .harness import MOCK_FILES, PART_SIZE, TRANSITION_WARNINGS, Grz
from .scenarios import SCENARIOS, Flow, Run, Scenario, expected_states

pytestmark = pytest.mark.e2e

_ARCHIVES = {"consented": services.ARCHIVE_CONSENTED, "non_consented": services.ARCHIVE_NON_CONSENTED}


@pytest.mark.parametrize("flow", ["process", "steps"])
@pytest.mark.parametrize("scenario", SCENARIOS, ids=lambda scenario: scenario.id)
def test_scenario(scenario: Scenario, flow: Flow, grz: Grz):
    for index, earlier in enumerate(scenario.earlier):
        grz.run_flow(flow, grz.upload(earlier, f"earlier-{index}"), (Run(),))
    grz.bfarm.received.clear()
    submission_id = grz.upload(scenario.submission, "submission")
    grz.run_flow(flow, submission_id, scenario.runs)

    outcome = scenario.outcome
    metadata = SubmissionMetadata(grz.work_dir / "submission" / "metadata" / "metadata.json")
    file_keys = {f"{submission_id}/files/{f.encrypted_file_path()}": f for f in metadata.files.values()}

    states, failure_reason = grz.states(submission_id)
    assert states == expected_states(scenario, flow)
    assert failure_reason == outcome.failure_reason

    for name, bucket in _ARCHIVES.items():
        archived = grz.keys(bucket, submission_id)
        if outcome.archive != name:
            assert not archived, f"{name} archive holds {archived}"
            continue
        assert set(file_keys) | {f"{submission_id}/metadata/metadata.json"} <= archived
        in_parts = False
        for key in file_keys:
            head = grz.s3().meta.client.head_object(Bucket=bucket, Key=key)
            in_parts |= "-" in head["ETag"]
            if head["ContentLength"] > PART_SIZE:
                assert "-" in head["ETag"], f"{key} is larger than a part but was not stored in parts"
        assert in_parts or not outcome.archived_in_parts, "no archived file was stored in parts"
        secret_key = crypt4gh.keys.get_private_key(MOCK_FILES / f"archive_{name}.sec", lambda: None)
        for key, file in file_keys.items():
            assert _decrypted_sha256(grz, bucket, key, secret_key) == file.file_checksum, key
    assert not grz.keys(services.ARCHIVE_INTERROGATION, submission_id), "interrogation bucket keeps staged copies"
    for bucket in services.BUCKETS:
        uploads = grz.s3().meta.client.list_multipart_uploads(Bucket=bucket, Prefix=f"{submission_id}/")
        assert not uploads.get("Uploads"), f"{bucket} keeps unfinished multipart uploads"

    inbox = grz.keys(services.INBOX, submission_id)
    if outcome.inbox_cleaned:
        assert inbox == {f"{submission_id}/cleaned", f"{submission_id}/metadata/metadata.json"}
    else:
        removed = {f"{submission_id}/files/{path}.c4gh" for path in scenario.submission.removed_from_inbox}
        assert set(file_keys) - removed <= inbox

    tan_g = metadata.content.submission.tan_g.lower()
    assert [r["SubmittedCase"]["tan"].lower() for r in grz.bfarm.received] == [tan_g] * outcome.pruefberichte

    submission = grz.submission(submission_id)
    if outcome.detailed_qc:
        assert submission.selected_for_qc is True
        assert len(grz.qc_output(flow, submission_id).read_text().splitlines()) == len(file_keys)
    else:
        assert not submission.selected_for_qc
        assert not (grz.qc_storage / submission_id / "files").exists(), "decrypted copies are left on local storage"

    warnings = [line for out in grz.outputs for line in out.splitlines() if any(w in line for w in TRANSITION_WARNINGS)]
    assert not warnings


def _decrypted_sha256(grz: Grz, bucket: str, key: str, secret_key: bytes) -> str:
    """Decrypt an archived object with the upstream crypt4gh library, not with grz-tools' own code."""
    body = io.BytesIO(grz.s3().Object(bucket, key).get()["Body"].read())
    plain = io.BytesIO()
    crypt4gh.lib.decrypt([(0, secret_key, None)], body, plain)
    return hashlib.sha256(plain.getvalue()).hexdigest()
