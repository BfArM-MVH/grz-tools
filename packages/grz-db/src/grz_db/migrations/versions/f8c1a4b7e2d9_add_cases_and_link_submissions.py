"""add cases table and link submissions

Revision ID: f8c1a4b7e2d9
Revises: c8d4a7f2b1e6
Create Date: 2026-07-08 00:00:00.000000+00:00

Introduces the ``cases`` table and links each submission to its case via ``submissions.case_id``.

A case's authoritative identity is ``psn`` (the RKI pseudonym), unique once assigned.
``submitter_id`` and ``local_case_id`` are resolution keys (indexed, but not a uniqueness
constraint) used to locate a case before a ``psn`` exists. Existing submissions are grouped by
``(submitter_id, pseudonym)`` (the ``pseudonym`` column holds the submitter-local case id) and
backfilled into cases. ``submissions.submitter_id`` is kept.

Also extends ``failurereasonenum`` with ``duplicate_initial``, recorded on a duplicate initial
submission: one being validated for a case whose initial submission already passed basic QC.
"""

from collections.abc import Sequence

import sqlalchemy as sa
from alembic import op
from sqlmodel.sql.sqltypes import AutoString

# revision identifiers, used by Alembic.
revision: str = "f8c1a4b7e2d9"
down_revision: str | Sequence[str] | None = "c8d4a7f2b1e6"
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None

# Mirrors grz_pydantic_models.submission.metadata.REDACTED_LOCAL_CASE_ID; hardcoded so the
# migration stays frozen. Redaction placeholders (this literal, or the empty string used by
# older archival code) are not case keys; SubmissionDb applies the same check at runtime.
PSEUDONYM_NON_KEYS = ["", "REDACTED_LOCAL_CASE_ID"]


def upgrade() -> None:
    """Upgrade schema."""
    bind = op.get_bind()
    submissions = sa.Table("submissions", sa.MetaData(), autoload_with=bind)

    # Fail fast only on genuinely unresolvable data: more than one QC-passed 'initial' submission
    # with the same (submitter_id, local_case_id). Pending or failed duplicates are legal, since a
    # case may have at most one initial submission that passed basic QC; they simply stay linked
    # to their case. Rows that cannot form a key are left unlinked.
    duplicate_initials = bind.execute(
        sa.select(
            submissions.c.submitter_id,
            submissions.c.pseudonym,
            sa.func.count().label("n"),
        )
        .where(
            submissions.c.pseudonym.is_not(None),
            submissions.c.pseudonym.not_in(PSEUDONYM_NON_KEYS),
            submissions.c.submitter_id.is_not(None),
            submissions.c.submission_type == "initial",
            submissions.c.basic_qc_passed.is_(True),
        )
        .group_by(submissions.c.submitter_id, submissions.c.pseudonym)
        .having(sa.func.count() > 1)
    ).fetchall()

    if duplicate_initials:
        groups = ", ".join(f"({row[0]}, {row[1]}): {row[2]}" for row in duplicate_initials)
        raise RuntimeError(
            "Cannot backfill cases: more than one QC-passed 'initial' submission shares the same "
            f"(submitter_id, local_case_id): {groups}"
        )

    # --- Create the cases table ---
    op.create_table(
        "cases",
        sa.Column("id", sa.Integer(), primary_key=True, autoincrement=True),
        sa.Column("psn", AutoString(), nullable=True),
        sa.Column("submitter_id", AutoString(), nullable=True),
        sa.Column("local_case_id", AutoString(), nullable=True),
    )
    # psn is the authoritative identity: unique when present.
    op.create_index(
        "ux_cases_psn",
        table_name="cases",
        columns=["psn"],
        unique=True,
        postgresql_where=sa.text("psn IS NOT NULL"),
        sqlite_where=sa.text("psn IS NOT NULL"),
    )
    # A submitter's local case ID identifies one patient, so it identifies one case. Unique only
    # where both halves are present: a case resolved by psn alone carries neither, and any number
    # of those may coexist. Without this two processes populating one patient's submissions both
    # find no case and both create one, which leaves the key ambiguous for good.
    op.create_index(
        "ux_cases_submitter_local_case",
        table_name="cases",
        columns=["submitter_id", "local_case_id"],
        unique=True,
        postgresql_where=sa.text("submitter_id IS NOT NULL AND local_case_id IS NOT NULL"),
        sqlite_where=sa.text("submitter_id IS NOT NULL AND local_case_id IS NOT NULL"),
    )

    # --- Link column on submissions ---
    op.add_column("submissions", sa.Column("case_id", sa.Integer(), nullable=True))
    op.create_index("ix_submissions_case_id", "submissions", ["case_id"])
    if bind.dialect.name == "postgresql":
        # SQLite cannot ALTER TABLE ADD CONSTRAINT; the ORM declares the FK regardless.
        op.create_foreign_key(
            "fk_submissions_case_id_cases",
            source_table="submissions",
            referent_table="cases",
            local_cols=["case_id"],
            remote_cols=["id"],
        )

    # --- Backfill: one case per distinct (submitter_id, pseudonym), storing the keys, then link ---
    # Lightweight, standalone Table objects (not reflected) describing just the columns we
    # need. A fresh MetaData per table avoids re-defining the ``submissions`` table reflected above.
    cases = sa.Table(
        "cases",
        sa.MetaData(),
        sa.Column("id", sa.Integer(), primary_key=True),
        sa.Column("submitter_id", AutoString()),
        sa.Column("local_case_id", AutoString()),
    )
    submissions_link = sa.Table(
        "submissions",
        sa.MetaData(),
        sa.Column("case_id", sa.Integer()),
        sa.Column("submitter_id", AutoString()),
        sa.Column("pseudonym", AutoString()),
        # declared as the existing native enum so comparisons render with the right type
        sa.Column(
            "submission_type",
            sa.Enum("initial", "followup", "addition", "correction", "test", name="submissiontype"),
        ),
    )

    # Only rows that can form a (submitter_id, pseudonym) key participate in the backfill.
    # Test submissions are never case-tracked; the explicit NULL check keeps rows with an unknown
    # type in the backfill, since SQL's three-valued logic would otherwise exclude them too.
    has_keys = sa.and_(
        submissions_link.c.pseudonym.is_not(None),
        submissions_link.c.pseudonym.not_in(PSEUDONYM_NON_KEYS),
        submissions_link.c.submitter_id.is_not(None),
        sa.or_(
            submissions_link.c.submission_type.is_(None),
            submissions_link.c.submission_type != "test",
        ),
    )
    # INSERT INTO cases (submitter_id, local_case_id)
    # SELECT DISTINCT submitter_id, pseudonym FROM submissions WHERE has_keys
    bind.execute(
        sa.insert(cases).from_select(
            ["submitter_id", "local_case_id"],
            sa.select(submissions_link.c.submitter_id, submissions_link.c.pseudonym).where(has_keys).distinct(),
        )
    )
    # Scalar subquery, one per submissions row, that looks up the matching case's id:
    #   SELECT c.id FROM cases c
    #   WHERE c.submitter_id = submissions.submitter_id AND c.local_case_id = submissions.pseudonym
    matching_case_id = (
        sa.select(cases.c.id)
        .where(
            cases.c.submitter_id == submissions_link.c.submitter_id,
            cases.c.local_case_id == submissions_link.c.pseudonym,
        )
        # mark `submissions_link` as coming from the enclosing UPDATE
        .correlate(submissions_link)
        # mark the SELECT as a single-value expression usable in `SET case_id = ...`.
        .scalar_subquery()
    )
    bind.execute(sa.update(submissions_link).values(case_id=matching_case_id).where(has_keys))

    # --- Enforce at most one QC-passed initial submission per case ---
    # Several initial submissions may be linked while pending basic QC (the data alone cannot
    # tell competing uploads apart); whichever passes basic QC first becomes the case's
    # QC-passed initial submission.
    op.create_index(
        "ux_submissions_one_initial_per_case",
        table_name="submissions",
        columns=["case_id"],
        unique=True,
        postgresql_where=sa.text("submission_type = 'initial' AND basic_qc_passed IS TRUE"),
        sqlite_where=sa.text("submission_type = 'initial' AND basic_qc_passed IS TRUE"),
    )

    # --- Failure reason for the duplicate initial submission ---
    # SQLite stores the enum as plain VARCHAR, so only PostgreSQL needs the type extended.
    # ADD VALUE inside a transaction requires PostgreSQL 12+, and even then the new value
    # cannot be used in the same transaction; this migration does not use it.
    if bind.dialect.name == "postgresql":
        op.execute("ALTER TYPE failurereasonenum ADD VALUE IF NOT EXISTS 'duplicate_initial'")


def downgrade() -> None:
    """Downgrade schema."""
    raise RuntimeError("Downgrades not supported.")
