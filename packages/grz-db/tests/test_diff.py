from grz_db.models.submission.diff import FieldDiff, SubmissionDiffCollection


def _collection() -> SubmissionDiffCollection:
    """A collection carrying one diff of every state, so filtering can be told apart from dropping."""
    collection = SubmissionDiffCollection()
    collection.append(FieldDiff.classify_field("submission_size", None, 42))  # added
    collection.append(FieldDiff.classify_field("submission_metadata", {"stale": True}, {"fresh": True}))  # updated
    collection.append(FieldDiff.classify_field("consented", True, False))  # updated
    collection.append(FieldDiff.classify_field("pseudonym", "old", None))  # deleted
    collection.append(FieldDiff.classify_field("tan_g", "same", "same"))  # unchanged
    return collection


def test_withhold_destructive_keeps_only_the_allowed_overwrites():
    committable, withheld = _collection().withhold_destructive({"submission_metadata"})

    assert [d.key for d in committable.updated] == ["submission_metadata"]
    assert [d.key for d in committable.deleted] == []
    assert sorted(d.key for d in withheld) == ["consented", "pseudonym"]


def test_withhold_destructive_always_keeps_additive_and_unchanged_diffs():
    """Writing a field that was NULL destroys nothing, so an allow-list must not hold it back."""
    committable, _ = _collection().withhold_destructive(set())

    assert [d.key for d in committable.added] == ["submission_size"]
    assert [d.key for d in committable.unchanged] == ["tan_g"]
    assert not committable.has_pending_destructive


def test_withhold_destructive_leaves_the_original_untouched():
    """The caller still reports on what was held back, so the source collection must not be mutated."""
    original = _collection()

    original.withhold_destructive({"submission_metadata"})

    assert len(list(original.pending)) == 4
    assert original.has_pending_destructive


def test_withhold_destructive_with_everything_allowed_changes_nothing():
    original = _collection()

    committable, withheld = original.withhold_destructive({d.key for d in original.pending})

    assert withheld == []
    assert [d.key for d in committable.pending] == [d.key for d in original.pending]
