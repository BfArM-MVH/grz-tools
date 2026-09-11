"""add download_error to failurereasonenum

Recorded when a submission cannot be downloaded from the inbox, the mirror image of the
existing ``upload_error``. Such failures were previously recorded as ``unknown``.

Revision ID: a7e3b1c94f20
Revises: c8d4a7f2b1e6
Create Date: 2026-08-19 00:00:00.000000+00:00
"""

from collections.abc import Sequence

from alembic import op

revision: str = "a7e3b1c94f20"
down_revision: str | Sequence[str] | None = "c8d4a7f2b1e6"
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None


def upgrade() -> None:
    """Upgrade schema."""
    # SQLite stores the enum as plain VARCHAR, so only PostgreSQL needs the type extended.
    # ADD VALUE inside a transaction requires PostgreSQL 12+, and even then the new value
    # cannot be used in the same transaction; this migration does not use it.
    if op.get_bind().dialect.name == "postgresql":
        op.execute("ALTER TYPE failurereasonenum ADD VALUE IF NOT EXISTS 'download_error'")


def downgrade() -> None:
    """Downgrade schema."""
    raise RuntimeError("Downgrades not supported.")
