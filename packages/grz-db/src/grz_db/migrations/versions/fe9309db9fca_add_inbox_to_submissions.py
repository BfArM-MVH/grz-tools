"""add inbox to submissions

Revision ID: fe9309db9fca
Revises: b7e4c1f9a2d3
Create Date: 2026-09-24 00:00:00.000000+00:00

Records which inbox a submission was downloaded from, so commands and the decryption key
lookup can find a submission again without being told the inbox again. The inbox's S3 bucket is
not stored: it follows from the configuration (``bucket or inbox_name``) at use time.
"""

from collections.abc import Sequence

import sqlalchemy as sa
from alembic import op

# revision identifiers, used by Alembic.
revision: str = "fe9309db9fca"
down_revision: str | Sequence[str] | None = "b7e4c1f9a2d3"
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None


def upgrade() -> None:
    """Upgrade schema."""
    op.add_column("submissions", sa.Column("inbox", sa.String(), nullable=True))


def downgrade() -> None:
    """Downgrade schema."""
    raise RuntimeError("Downgrades not supported.")
