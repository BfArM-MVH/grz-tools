from grz_db.models.submission import SubmissionType
from grz_db.models.submission.diff import (
    CaseLinkDiff,
    DiffState,
    DonorDiff,
    DonorsDiffCollection,
    FieldDiff,
    SubmissionChangeSet,
    SubmissionDiffCollection,
)


def _collection() -> SubmissionDiffCollection:
    """A collection carrying one diff of every state, so filtering can be told apart from dropping."""
    collection = SubmissionDiffCollection()
    collection.append(FieldDiff.classify_field("submission_size", None, 42))  # added
    collection.append(FieldDiff.classify_field("submission_metadata", {"stale": True}, {"fresh": True}))  # updated
    collection.append(FieldDiff.classify_field("consented", True, False))  # updated
    collection.append(FieldDiff.classify_field("local_case_id", "old", None))  # deleted
    collection.append(FieldDiff.classify_field("tan_g", "same", "same"))  # unchanged
    return collection


def test_withhold_destructive_keeps_only_the_allowed_overwrites():
    committable, withheld = _collection().withhold_destructive({"submission_metadata"})

    assert [d.key for d in committable.updated] == ["submission_metadata"]
    assert [d.key for d in committable.deleted] == []
    assert sorted(d.key for d in withheld) == ["consented", "local_case_id"]


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


def _donor(pseudonym: str, state: DiffState) -> DonorDiff:
    return DonorDiff(before=None, after=None, state=state, pseudonym=pseudonym)


def _case_link(before: int | None) -> CaseLinkDiff:
    return CaseLinkDiff(
        before=before, after=2, submitter_id="123456789", local_case_id="case", submission_type=SubmissionType.initial
    )


def _change_set(case_link: CaseLinkDiff | None = None) -> SubmissionChangeSet:
    """A change set carrying every kind of change, so what stays can be told apart from what is held back."""
    donors = DonorsDiffCollection()
    donors.append(_donor("index", DiffState.NEW))
    donors.append(_donor("mother", DiffState.UPDATED))
    donors.append(_donor("father", DiffState.DELETED))
    donors.append(_donor("sister", DiffState.UNCHANGED))
    return SubmissionChangeSet(fields=_collection(), donors=donors, case_link=case_link)


def test_change_set_withholds_updated_and_deleted_donors_but_adds_new_ones():
    committable, withheld = _change_set().withhold_destructive(set())

    assert [d.pseudonym for d in committable.donors.pending] == ["index"]
    assert [d.pseudonym for d in committable.donors.unchanged] == ["sister"]
    assert [d.pseudonym for d in withheld.donors.pending] == ["mother", "father"]


def test_change_set_writes_donor_changes_that_allowed_names():
    committable, withheld = _change_set().withhold_destructive({"donors"})

    assert [d.pseudonym for d in committable.donors.pending] == ["index", "mother", "father"]
    assert not withheld.donors.has_pending


def test_change_set_withholds_fields_like_the_field_collection():
    committable, withheld = _change_set().withhold_destructive({"submission_metadata"})

    assert [d.key for d in committable.fields.pending] == ["submission_size", "submission_metadata"]
    assert sorted(d.key for d in withheld.fields.pending) == ["consented", "local_case_id"]


def test_change_set_withholds_a_changed_case_link_unless_allowed():
    link = _case_link(before=1)

    committable, withheld = _change_set(link).withhold_destructive(set())
    assert committable.case_link is None
    assert withheld.case_link is link

    committable, withheld = _change_set(link).withhold_destructive({"case_id"})
    assert committable.case_link is link
    assert withheld.case_link is None


def test_change_set_always_keeps_a_first_case_link():
    """Linking a submission that has no case yet destroys nothing."""
    link = _case_link(before=None)

    committable, withheld = _change_set(link).withhold_destructive(set())

    assert committable.case_link is link
    assert withheld.case_link is None


def test_change_set_withheld_part_names_exactly_what_was_held_back():
    original = _change_set(_case_link(before=1))

    committable, withheld = original.withhold_destructive(set())

    assert not committable.has_pending_destructive
    assert sorted(withheld.destructive_changes) == sorted(original.destructive_changes)


def test_change_set_withholding_leaves_the_original_untouched():
    original = _change_set(_case_link(before=1))
    destructive_before = original.destructive_changes

    original.withhold_destructive(set())

    assert original.destructive_changes == destructive_before
