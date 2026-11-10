"""add 'processing' and 'processed' states

Revision ID: 2f6a24f8db05
Revises: f8c1a4b7e2d9
Create Date: 2026-09-17 00:00:00.000000+00:00

Adds the states that ``grzctl process`` records to the PostgreSQL ``submissionstateenum`` type.
"""

from collections.abc import Sequence

from alembic import op

revision: str = "2f6a24f8db05"
down_revision: str | Sequence[str] | None = "f8c1a4b7e2d9"
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None


def upgrade() -> None:
    bind = op.get_bind()
    # SQLite stores the enum as plain VARCHAR, so only PostgreSQL needs the type extended.
    # ADD VALUE inside a transaction requires PostgreSQL 12+, and even then the new values
    # cannot be used in the same transaction; this migration does not use them.
    if bind.dialect.name == "postgresql":
        op.execute("ALTER TYPE submissionstateenum ADD VALUE IF NOT EXISTS 'PROCESSING' BEFORE 'UPLOADING'")
        op.execute("ALTER TYPE submissionstateenum ADD VALUE IF NOT EXISTS 'PROCESSED' BEFORE 'FINISHED'")


def downgrade() -> None:
    """Downgrade schema."""
    raise RuntimeError("Downgrades not supported.")
