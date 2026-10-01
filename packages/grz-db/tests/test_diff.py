from grz_db.models.submission import SubmissionType
from grz_db.models.submission.diff import (
    CASE_LINK_KEY,
    DONORS_KEY,
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


def _donor(pseudonym: str, state: DiffState) -> DonorDiff:
    return DonorDiff(before=None, after=None, state=state, pseudonym=pseudonym)


def _case_link(before: int | None) -> CaseLinkDiff:
    return CaseLinkDiff(
        before=before, after=2, submitter_id="123456789", local_case_id="case", submission_type=SubmissionType.initial
    )


def _change_set(case_link: CaseLinkDiff | None = None) -> SubmissionChangeSet:
    """A change set carrying every kind of change, so what counts can be told apart from what does not."""
    donors = DonorsDiffCollection()
    donors.append(_donor("index", DiffState.NEW))
    donors.append(_donor("mother", DiffState.UPDATED))
    donors.append(_donor("father", DiffState.DELETED))
    donors.append(_donor("sister", DiffState.UNCHANGED))
    return SubmissionChangeSet(fields=_collection(), donors=donors, case_link=case_link)


def test_every_destructive_change_is_named_when_nothing_is_allowed():
    undeclared = _change_set(_case_link(before=1)).undeclared_destructive_changes()

    assert undeclared == [
        "case link (case 1)",
        "consented",
        "donor 'father'",
        "donor 'mother'",
        "local_case_id",
        "submission_metadata",
    ]


def test_additive_changes_are_never_named():
    """Writing a field that was NULL, or a donor the database does not have, destroys nothing."""
    undeclared = _change_set().undeclared_destructive_changes()

    assert "submission_size" not in undeclared
    assert "donor 'index'" not in undeclared


def test_a_field_the_allow_list_names_drops_out():
    undeclared = _change_set().undeclared_destructive_changes({"submission_metadata"})

    assert undeclared == ["consented", "donor 'father'", "donor 'mother'", "local_case_id"]


def test_the_donors_key_covers_every_donor_update_and_delete():
    undeclared = _change_set().undeclared_destructive_changes({DONORS_KEY})

    assert undeclared == ["consented", "local_case_id", "submission_metadata"]


def test_the_case_link_key_covers_a_changed_link():
    undeclared = _change_set(_case_link(before=1)).undeclared_destructive_changes({CASE_LINK_KEY})

    assert "case link (case 1)" not in undeclared


def test_a_first_case_link_is_not_destructive():
    """Linking a submission that has no case yet destroys nothing."""
    undeclared = _change_set(_case_link(before=None)).undeclared_destructive_changes()

    assert not [name for name in undeclared if name.startswith("case link")]


def test_allowing_everything_leaves_nothing_undeclared():
    change_set = _change_set(_case_link(before=1))
    allowed = {d.key for d in change_set.fields.pending} | {DONORS_KEY, CASE_LINK_KEY}

    assert change_set.undeclared_destructive_changes(allowed) == []


def test_destructive_changes_names_what_an_empty_allow_list_leaves():
    change_set = _change_set(_case_link(before=1))

    assert change_set.destructive_changes == change_set.undeclared_destructive_changes()
    assert change_set.has_pending_destructive


def test_asking_leaves_the_change_set_untouched():
    """The caller reports on what it refused, so the change set it reports from must survive."""
    change_set = _change_set(_case_link(before=1))
    pending_before = list(change_set.pending_changes)

    change_set.undeclared_destructive_changes({DONORS_KEY})

    assert change_set.pending_changes == pending_before
    assert change_set.destructive_changes
