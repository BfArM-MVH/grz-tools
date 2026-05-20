from __future__ import annotations

from collections.abc import Generator
from dataclasses import dataclass, field
from enum import StrEnum
from typing import TYPE_CHECKING, Self, cast

if TYPE_CHECKING:
    from grz_db.models.submission import Donor as DbDonor


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
class SubmissionDiffCollection:
    """Holds the result of diffing submission-level metadata against the database.

    Fields are categorised the same way as :class:`DonorDiff`:

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
        """All diffs that need to be written to the database (added + updated + deleted)."""
        yield from self.added
        yield from self.updated
        yield from self.deleted

    @property
    def has_pending(self) -> bool:
        """True if any field needs to be written to the database (added, updated, or deleted)."""
        return len(self.added) > 0 or len(self.updated) > 0 or len(self.deleted) > 0

    @property
    def has_pending_destructive(self) -> bool:
        """True if any field will overwrite or remove an existing database value (updated or deleted)."""
        return len(self.updated) > 0 or len(self.deleted) > 0

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
                donor.changes.append(
                    FieldDiff.classify_field(str(f), getattr(old_value, str(f), None), getattr(new_value, str(f), None))
                )

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
        """All diffs that need to be written to the database (added + updated + deleted)."""
        yield from self.added
        yield from self.updated
        yield from self.deleted

    @property
    def has_pending(self) -> bool:
        """True if any donor needs to be written to the database (added, updated, or deleted)."""
        return len(self.added) > 0 or len(self.updated) > 0 or len(self.deleted) > 0

    @property
    def has_pending_destructive(self) -> bool:
        """True if any donor will overwrite or remove an existing database record (updated or deleted)."""
        return len(self.updated) > 0 or len(self.deleted) > 0

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
