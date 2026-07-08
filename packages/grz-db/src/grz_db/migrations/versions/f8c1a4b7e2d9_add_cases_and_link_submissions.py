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

    # Fail fast only on genuinely inconsistent data: more than one 'initial' sharing a
    # (submitter_id, local_case_id). Rows that cannot form a key are simply left unlinked.
    duplicate_initials = bind.execute(
        sa.text(
            "SELECT submitter_id, pseudonym, COUNT(*) AS n FROM submissions "
            "WHERE pseudonym IS NOT NULL AND submitter_id IS NOT NULL AND submission_type = 'initial' "
            "GROUP BY submitter_id, pseudonym HAVING COUNT(*) > 1"
        )
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
        "cases",
        ["psn"],
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
            "submissions",
            "cases",
            ["case_id"],
            ["id"],
        )

    # --- Backfill: one case per distinct (submitter_id, pseudonym), storing the keys, then link ---
    op.execute(
        "INSERT INTO cases (submitter_id, local_case_id) "
        "SELECT DISTINCT submitter_id, pseudonym FROM submissions "
        "WHERE pseudonym IS NOT NULL AND submitter_id IS NOT NULL"
    )
    op.execute(
        "UPDATE submissions SET case_id = ("
        "SELECT c.id FROM cases c "
        "WHERE c.submitter_id = submissions.submitter_id AND c.local_case_id = submissions.pseudonym"
        ") WHERE pseudonym IS NOT NULL AND submitter_id IS NOT NULL"
    )

    # --- Enforce at most one initial submission per case ---
    op.create_index(
        "ux_submissions_one_initial_per_case",
        "submissions",
        ["case_id"],
        unique=True,
        postgresql_where=sa.text("submission_type = 'initial'"),
        sqlite_where=sa.text("submission_type = 'initial'"),
    )


def downgrade() -> None:
    """Downgrade schema."""
    raise RuntimeError("Downgrades not supported.")
