from __future__ import annotations

from collections.abc import Container, Generator
from dataclasses import dataclass, field
from enum import StrEnum
from typing import TYPE_CHECKING, Self, cast

if TYPE_CHECKING:
    from grz_pydantic_models.submission.metadata import SubmissionType

    from grz_db.errors import AmbiguousCaseError
    from grz_db.models.submission import Donor as DbDonor


DONORS_KEY = "donors"
"""Names every donor update and delete at once in an allow-list.

A donor row is keyed by the pseudonym a change may itself be replacing, so there is no finer
unit to name.
"""

CASE_LINK_KEY = "case_id"
"""Names a changed case link, both in an allow-list and in ``ignore_fields``.

It is the column the link is stored in, so operators name it the way they name a field.
"""


class DiffState(StrEnum):
    NEW = "new"
    UPDATED = "update"
    DELETED = "delete"
    UNCHANGED = "unchanged"


@dataclass
class Diff[T]:
    """A single difference between old and new value.

    :param before: Old value
    :param after: New value
    :param state: State of the change
    """

    before: T | None
    after: T | None
    state: DiffState

    @classmethod
    def classify(
        cls,
        old_value: T | None,
        new_value: T | None,
    ) -> Self:
        """Classify the difference between old and new value.

        :param old_value: Old value (``None`` if unset).
        :param new_value: New value (``None`` if unset).
        """
        match (old_value, new_value):
            case (None, None):
                state = DiffState.UNCHANGED
            case (None, _):
                state = DiffState.NEW
            case (_, None):
                state = DiffState.DELETED
            case _ if old_value != new_value:
                state = DiffState.UPDATED
            case _:
                state = DiffState.UNCHANGED
        return cls(old_value, new_value, state)


@dataclass
class FieldDiff[T]:
    """A single field-level difference between old and new value.

    :param key: Name of the field that changed.
    :param diff: Difference between old and new value, including the state of the change.
    """

    key: str
    diff: Diff[T]

    @classmethod
    def classify_field(
        cls,
        key: str,
        old_value: T | None,
        new_value: T | None,
    ) -> Self:
        """Classify the difference between old and new value.

        :param key: Field name being compared.
        :param old_value: Old value (``None`` if unset).
        :param new_value: New value (``None`` if unset).
        """
        return cls(key, Diff.classify(old_value, new_value))


@dataclass
class CaseLinkDiff:
    """A pending change to a submission's case link, resolved during :meth:`SubmissionDb.diff`.

    Only constructed when committing would change the link, so its mere presence on a
    :class:`SubmissionChangeSet` means there is something to write.

    :param before: ``case_id`` currently stored on the submission, if any.
    :param after: Primary key of the matching existing case, or ``None`` if committing
        will create a new case.
    :param submitter_id: Submitter identifier used to resolve the case.
    :param local_case_id: Submitter-local case identifier used to resolve the case.
    :param submission_type: Type of the submission. Checked for case-trackability (``test`` is
        rejected), not against the case.
    :param psn: RKI pseudonym used to resolve the case, when the caller knows one. Carried so
        that :meth:`SubmissionDb.commit_changes` can store it on a case it creates, and so a
        ``PsnResolver`` reaches the same value the diff resolved with. Defaults to ``None``,
        since no psn is derivable from a submitter's metadata.
    """

    before: int | None
    after: int | None
    submitter_id: str | None
    local_case_id: str | None
    submission_type: SubmissionType
    psn: str | None = None

    @property
    def state(self) -> DiffState:
        """NEW when the submission is not yet linked, UPDATED when an existing link would change.

        DELETED cannot occur: a diff never proposes unlinking; ``after is None`` means a new
        case would be created instead.
        """
        return DiffState.NEW if self.before is None else DiffState.UPDATED


@dataclass
class SubmissionDiffCollection:
    """Holds the result of diffing submission-level metadata against the database.

    Fields are categorised the same way as :class:`DonorsDiffCollection`:

    :param added: Fields that were ``None`` in the database and now have a value.
    :param updated: Fields whose non-null database value differs from the new value.
    :param deleted: Fields that had a value in the database but are now ``None``.
    :param unchanged: Fields whose value is already in sync with the new value.
    """

    added: list[FieldDiff] = field(default_factory=list)
    updated: list[FieldDiff] = field(default_factory=list)
    deleted: list[FieldDiff] = field(default_factory=list)
    unchanged: list[FieldDiff] = field(default_factory=list)

    @property
    def pending(self) -> Generator[FieldDiff, None, None]:
        """All field diffs that need to be written to the database (added + updated + deleted)."""
        yield from self.added
        yield from self.updated
        yield from self.deleted

    @property
    def has_pending(self) -> bool:
        """True if any field needs to be written to the database (added, updated, or deleted)."""
        return len(self.added) > 0 or len(self.updated) > 0 or len(self.deleted) > 0

    def append(self, field_diff: FieldDiff):
        match field_diff.diff.state:
            case DiffState.UPDATED:
                self.updated.append(field_diff)
            case DiffState.NEW:
                self.added.append(field_diff)
            case DiffState.DELETED:
                self.deleted.append(field_diff)
            case DiffState.UNCHANGED:
                self.unchanged.append(field_diff)


@dataclass
class DonorDiff(Diff["DbDonor"]):
    """Holds the diff result for a single donor.

    :param diff: The new/updated :class:`Donor` model to be written.
    :param changes: List of :class:`FieldDiff` instances for every field that differs.
    :param state: Classification of the donor relative to the current database state.
    """

    pseudonym: str | None = None
    changes: list[FieldDiff] = field(default_factory=list)

    @classmethod
    def classify(cls, old_value: DbDonor | None, new_value: DbDonor | None) -> Self:
        # lazy import to avoid import cycle
        from grz_db.models.submission import Donor as DbDonor  # noqa: PLC0415

        # Keep `result` as `Self` so the return type is correct; use `donor` (cast to the
        # concrete DonorDiff) for attribute access that the parent class doesn't expose.
        result = super().classify(old_value, new_value)
        donor = cast("DonorDiff", result)

        if donor.state != DiffState.UNCHANGED:
            for f in sorted(DbDonor.model_fields.keys() - {"submission_id", "pseudonym"}):
                field_diff = FieldDiff.classify_field(
                    str(f), getattr(old_value, str(f), None), getattr(new_value, str(f), None)
                )
                if field_diff.diff.state != DiffState.UNCHANGED:
                    donor.changes.append(field_diff)

        match (old_value, new_value):
            case (None, None):
                donor.pseudonym = None
            case (None, _):
                donor.pseudonym = cast("DbDonor", new_value).pseudonym
            case (_, None):
                donor.pseudonym = cast("DbDonor", old_value).pseudonym
            case (_, _):
                old_donor = cast("DbDonor", old_value)
                new_donor = cast("DbDonor", new_value)
                if old_donor.pseudonym != new_donor.pseudonym:
                    raise ValueError(f"Pseudonym mismatch: '{old_donor.pseudonym}' != '{new_donor.pseudonym}'")
                donor.pseudonym = old_donor.pseudonym

        return result


@dataclass
class DonorsDiffCollection:
    """Aggregated result of comparing all donors in a metadata file to the database.

    :param added: Donors present in metadata but not yet in the database.
    :param updated: Donors whose metadata differs from the database record.
    :param deleted: Donors present in the database but absent from the metadata.
    :param unchanged: Donors whose metadata is already in sync with the database.
    """

    added: list[DonorDiff] = field(default_factory=list)
    updated: list[DonorDiff] = field(default_factory=list)
    deleted: list[DonorDiff] = field(default_factory=list)
    unchanged: list[DonorDiff] = field(default_factory=list)

    @property
    def pending(self) -> Generator[DonorDiff, None, None]:
        """All donor diffs that need to be written to the database (added + updated + deleted)."""
        yield from self.added
        yield from self.updated
        yield from self.deleted

    @property
    def has_pending(self) -> bool:
        """True if any donor needs to be written to the database (added, updated, or deleted)."""
        return len(self.added) > 0 or len(self.updated) > 0 or len(self.deleted) > 0

    def append(self, donor_diff: DonorDiff):
        match donor_diff.state:
            case DiffState.UPDATED:
                self.updated.append(donor_diff)
            case DiffState.NEW:
                self.added.append(donor_diff)
            case DiffState.DELETED:
                self.deleted.append(donor_diff)
            case DiffState.UNCHANGED:
                self.unchanged.append(donor_diff)


@dataclass
class SubmissionChangeSet:
    """Everything committing a metadata file would change for one submission.

    :param fields: Column-level diffs on the submission row.
    :param donors: Donor rows to add, update, or delete.
    :param case_link: Pending change to the submission's case link, or ``None`` if the
        link is already in sync (or case resolution was skipped). Applied via
        :meth:`SubmissionDb.assign_case`, not as a column write.
    :param case_link_error: Why the case could not be resolved, when it could not be.
        Whether a case can be identified says nothing about whether the rest of the
        metadata should be recorded, so this is reported alongside the other diffs
        rather than raised: the caller decides what an unlinkable submission is worth.
        ``case_link`` is ``None`` whenever this is set, since there is nothing to write.
    """

    fields: SubmissionDiffCollection = field(default_factory=SubmissionDiffCollection)
    donors: DonorsDiffCollection = field(default_factory=DonorsDiffCollection)
    case_link: CaseLinkDiff | None = None
    case_link_error: AmbiguousCaseError | None = None

    @property
    def has_pending(self) -> bool:
        """True if committing would write anything to the database."""
        return self.fields.has_pending or self.donors.has_pending or self.case_link is not None

    @property
    def pending_changes(self) -> list[str]:
        """Names of every pending change: field diffs, donor diffs, and the case link."""
        changes = [d.key for d in self.fields.pending]
        changes += [f"donor '{d.pseudonym}'" for d in self.donors.pending]
        if self.case_link is not None:
            source = f"case {self.case_link.before}" if self.case_link.before is not None else "unlinked"
            target = f"case {self.case_link.after}" if self.case_link.after is not None else "new case"
            changes.append(f"case link ({source} -> {target})")
        return sorted(changes)

    @property
    def destructive_changes(self) -> list[str]:
        """Names of the existing database values committing would overwrite or remove.

        An existing case link counts: reverting it would undo a deliberate ``case relink``.
        """
        return self.undeclared_destructive_changes()

    def undeclared_destructive_changes(self, allowed: Container[str] = frozenset()) -> list[str]:
        """Names of the destructive changes ``allowed`` does not cover.

        Additive changes never appear: a filled field, a new donor and a first case link destroy
        nothing. An empty list means the whole change set may be committed as it stands.

        :param allowed: What may be overwritten or removed. A field matches by key, every donor
            update and delete together as :data:`DONORS_KEY`, and a changed case link as
            :data:`CASE_LINK_KEY`.
        :returns: One name per change, sorted.
        """
        changes = [d.key for d in (*self.fields.updated, *self.fields.deleted) if d.key not in allowed]
        if DONORS_KEY not in allowed:
            changes += [f"donor '{d.pseudonym}'" for d in (*self.donors.updated, *self.donors.deleted)]
        if self.case_link is not None and self.case_link.state is DiffState.UPDATED and CASE_LINK_KEY not in allowed:
            changes.append(f"case link (case {self.case_link.before})")
        return sorted(changes)

    @property
    def has_pending_destructive(self) -> bool:
        """True if committing would overwrite or remove any existing database value."""
        return bool(self.destructive_changes)
