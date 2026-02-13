"""modify submissions with new column size and metadata

Revision ID: 834bd50b8734
Revises: 8aa6fc8d118a
Create Date: 2026-02-12 14:58:52.096841+00:00

"""

from collections.abc import Sequence

import sqlalchemy as sa
from alembic import op


# revision identifiers, used by Alembic.
revision: str = '834bd50b8734'
down_revision: str | Sequence[str] | None = '8aa6fc8d118a'
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None


def upgrade() -> None:
    """Upgrade schema."""
    op.add_column("submissions", sa.Column("submission_size", sa.Numeric(10, 2), nullable=True))
    op.add_column("submissions", sa.Column("submission_metadata", sa.Text(), nullable=True))



def downgrade() -> None:
    """Downgrade schema."""
    raise RuntimeError("Downgrades not supported.")
