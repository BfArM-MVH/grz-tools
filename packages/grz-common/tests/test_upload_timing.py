"""Tests for deriving the upload timing of a submission from its inbox objects."""

import datetime
from types import SimpleNamespace
from unittest.mock import patch

from grz_common.transfer import upload_timing
from grz_common.workers.download import query_submissions

SUBMISSION_ID = "260914050_2025-09-15_c64603a7"
T0 = datetime.datetime(2025, 9, 15, 8, 0, tzinfo=datetime.UTC)


def _at(minutes: int) -> datetime.datetime:
    return T0 + datetime.timedelta(minutes=minutes)


def test_upload_timing_spans_marker_to_last_object():
    """The ``version`` marker starts the upload, the newest object ends it."""
    timing = upload_timing(
        SUBMISSION_ID,
        {
            f"{SUBMISSION_ID}/version": _at(0),
            f"{SUBMISSION_ID}/files/a.c4gh": _at(10),
            f"{SUBMISSION_ID}/files/b.c4gh": _at(25),
            f"{SUBMISSION_ID}/metadata/metadata.json": _at(26),
        },
    )
    assert timing == (_at(0), _at(26), False)
    assert timing.duration == datetime.timedelta(minutes=26)


def test_upload_timing_detects_a_resumed_upload():
    """A resumed upload rewrites the marker, so an object older than the marker gives it away."""
    timing = upload_timing(
        SUBMISSION_ID,
        {
            f"{SUBMISSION_ID}/files/a.c4gh": _at(10),  # from the interrupted first run
            f"{SUBMISSION_ID}/version": _at(60),
            f"{SUBMISSION_ID}/files/b.c4gh": _at(70),
            f"{SUBMISSION_ID}/metadata/metadata.json": _at(71),
        },
    )
    assert timing == (_at(60), _at(71), True)


def test_upload_timing_without_marker_cannot_tell():
    """Without the marker (older grz-cli), the oldest object starts the upload and resumption is unknown."""
    timing = upload_timing(
        SUBMISSION_ID,
        {f"{SUBMISSION_ID}/files/a.c4gh": _at(10), f"{SUBMISSION_ID}/metadata/metadata.json": _at(26)},
    )
    assert timing == (_at(10), _at(26), None)


def test_upload_timing_needs_an_object():
    """An inbox without uploaded objects (cleaned, or only the marker) has no timing; it is not the LE's fault."""
    assert upload_timing(SUBMISSION_ID, {f"{SUBMISSION_ID}/version": _at(0)}) is None
    assert upload_timing(SUBMISSION_ID, {}) is None


def test_query_submissions_takes_the_duration_from_the_upload_timing():
    """``grzctl list`` shows the same duration as the monitoring report, not the oldest to the newest object."""
    only_marker = "260914050_2025-09-16_00000000"
    objects = [
        {"Key": f"{SUBMISSION_ID}/files/a.c4gh", "LastModified": _at(10), "Size": 1},
        {"Key": f"{SUBMISSION_ID}/version", "LastModified": _at(60), "Size": 1},
        {"Key": f"{SUBMISSION_ID}/metadata/metadata.json", "LastModified": _at(71), "Size": 1},
        {"Key": f"{only_marker}/version", "LastModified": _at(80), "Size": 1},
    ]
    with patch("grz_common.workers.download.init_s3_client") as init_s3_client:
        init_s3_client.return_value.get_paginator.return_value.paginate.return_value = [{"Contents": objects}]
        summaries = query_submissions(SimpleNamespace(bucket="inbox"), show_cleaned=False)  # type: ignore[arg-type]

    durations = {summary.submission_id: summary.upload_duration for summary in summaries}
    assert durations == {SUBMISSION_ID: datetime.timedelta(minutes=11), only_marker: None}
