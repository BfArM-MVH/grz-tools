"""add 'deleting' and 'deleted' states

Revision ID: b70d8efaf2bb
Revises: f49beda8afac
Create Date: 2026-11-19 00:00:00.000000+00:00

Adds the states that ``grzctl delete`` records to the PostgreSQL ``submissionstateenum`` type.
"""

from collections.abc import Sequence

from alembic import op

revision: str = "b70d8efaf2bb"
down_revision: str | Sequence[str] | None = "f49beda8afac"
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None


def upgrade() -> None:
    bind = op.get_bind()
    # SQLite stores the enum as plain VARCHAR, so only PostgreSQL needs the type extended.
    # ADD VALUE inside a transaction requires PostgreSQL 12+, and even then the new values
    # cannot be used in the same transaction; this migration does not use them.
    if bind.dialect.name == "postgresql":
        op.execute("ALTER TYPE submissionstateenum ADD VALUE IF NOT EXISTS 'DELETING' BEFORE 'ERROR'")
        op.execute("ALTER TYPE submissionstateenum ADD VALUE IF NOT EXISTS 'DELETED' BEFORE 'ERROR'")


def downgrade() -> None:
    """Downgrade schema."""
    raise RuntimeError("Downgrades not supported.")
