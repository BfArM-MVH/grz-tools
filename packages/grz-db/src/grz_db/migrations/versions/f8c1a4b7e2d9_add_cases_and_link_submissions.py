"""add cases table and link submissions

Revision ID: f8c1a4b7e2d9
Revises: 66f36abbea34
Create Date: 2026-07-08 00:00:00.000000+00:00

Introduces the ``cases`` table and links each submission to its case via ``submissions.case_id``.

A case's authoritative identity is ``psn`` (the RKI pseudonym), unique once assigned.
``submitter_id`` and ``local_case_id`` are resolution keys (indexed, but not a uniqueness
constraint) used to locate a case before a ``psn`` exists. Existing submissions are grouped by
``(submitter_id, pseudonym)`` and backfilled into cases. ``submissions.submitter_id`` is kept.
"""

from collections.abc import Sequence

import sqlalchemy as sa
from alembic import op
from sqlmodel.sql.sqltypes import AutoString

# revision identifiers, used by Alembic.
revision: str = "f8c1a4b7e2d9"
down_revision: str | Sequence[str] | None = "66f36abbea34"
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None


def upgrade() -> None:
    """Upgrade schema."""
    bind = op.get_bind()
    submissions = sa.Table("submissions", sa.MetaData(), autoload_with=bind)

    # Fail fast only on genuinely inconsistent data: more than one 'initial' sharing a
    # (submitter_id, local_case_id). Rows that cannot form a key are simply left unlinked.
    duplicate_initials = bind.execute(
        sa.select(
            submissions.c.submitter_id,
            submissions.c.pseudonym,
            sa.func.count().label("n"),
        )
        .where(
            submissions.c.pseudonym.is_not(None),
            submissions.c.submitter_id.is_not(None),
            submissions.c.submission_type == "initial",
        )
        .group_by(submissions.c.submitter_id, submissions.c.pseudonym)
        .having(sa.func.count() > 1)
    ).fetchall()

    if duplicate_initials:
        groups = ", ".join(f"({row[0]}, {row[1]}): {row[2]}" for row in duplicate_initials)
        raise RuntimeError(
            "Cannot backfill cases: more than one 'initial' submission shares the same "
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
    # (submitter_id, local_case_id) is a lookup key, not a uniqueness constraint.
    op.create_index("ix_cases_submitter_local_case", "cases", ["submitter_id", "local_case_id"])

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
    # need. Using a separate MetaData per table avoids clashing with the real ORM metadata.
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
    )

    # Only rows that can form a (submitter_id, pseudonym) key participate in the backfill.
    has_keys = sa.and_(
        submissions_link.c.pseudonym.is_not(None),
        submissions_link.c.submitter_id.is_not(None),
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
    # UPDATE submissions SET case_id = (matching_case_id) WHERE has_keys
    bind.execute(sa.update(submissions_link).values(case_id=matching_case_id).where(has_keys))

    # --- Enforce at most one initial submission per case ---
    op.create_index(
        "ux_submissions_one_initial_per_case",
        table_name="submissions",
        columns=["case_id"],
        unique=True,
        postgresql_where=sa.text("submission_type = 'initial'"),
        sqlite_where=sa.text("submission_type = 'initial'"),
    )


def downgrade() -> None:
    """Downgrade schema."""
    raise RuntimeError("Downgrades not supported.")
