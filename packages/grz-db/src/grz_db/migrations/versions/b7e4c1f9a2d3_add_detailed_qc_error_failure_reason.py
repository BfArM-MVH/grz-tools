"""add 'detailed_qc_error' failure reason

Revision ID: b7e4c1f9a2d3
Revises: 2f6a24f8db05
Create Date: 2026-09-17 00:00:00.000000+00:00

Adds the reason that ``grzctl process`` records when the configured detailed QC workflow fails,
to the PostgreSQL ``failurereasonenum`` type.
"""

from collections.abc import Sequence

from alembic import op

revision: str = "b7e4c1f9a2d3"
down_revision: str | Sequence[str] | None = "2f6a24f8db05"
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None


def upgrade() -> None:
    bind = op.get_bind()
    # SQLite stores the enum as plain VARCHAR, so only PostgreSQL needs the type extended.
    # ADD VALUE inside a transaction requires PostgreSQL 12+, and even then the new value
    # cannot be used in the same transaction; this migration does not use it.
    if bind.dialect.name == "postgresql":
        op.execute("ALTER TYPE failurereasonenum ADD VALUE IF NOT EXISTS 'detailed_qc_error'")


def downgrade() -> None:
    """Downgrade schema."""
    raise RuntimeError("Downgrades not supported.")
