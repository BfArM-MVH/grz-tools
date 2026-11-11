"""add upload timing to submissions

Revision ID: a7c3e9d1f5b2
Revises: fe9309db9fca
Create Date: 2026-10-29 00:00:00.000000+00:00

Records when the inbox received the first and the last object of a submission's upload, and
whether the upload was resumed after an interruption. ``grzctl download --populate`` fills the
columns from the inbox's object listing, which ``grzctl clean`` erases, so rows of submissions
cleaned before this migration stay NULL. ``grzctl report monitoring`` aggregates transfer speeds
per Leistungserbringer and month from them.
"""

from collections.abc import Sequence

import sqlalchemy as sa
from alembic import op

# revision identifiers, used by Alembic.
revision: str = "a7c3e9d1f5b2"
down_revision: str | Sequence[str] | None = "fe9309db9fca"
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None


def upgrade() -> None:
    """Upgrade schema."""
    op.add_column("submissions", sa.Column("upload_started_at", sa.DateTime(timezone=True), nullable=True))
    op.add_column("submissions", sa.Column("upload_finished_at", sa.DateTime(timezone=True), nullable=True))
    op.add_column("submissions", sa.Column("upload_resumed", sa.Boolean(), nullable=True))


def downgrade() -> None:
    """Downgrade schema."""
    raise RuntimeError("Downgrades not supported.")
