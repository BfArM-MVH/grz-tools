"""add failure_reason to submission_states

Revision ID: 524b2052f143
Revises: f3a9c7d2e481
Create Date: 2026-04-20 10:32:47.329276+00:00

"""

from collections.abc import Sequence

import sqlalchemy as sa
from alembic import op

# revision identifiers, used by Alembic.
revision: str = "524b2052f143"
down_revision: str | Sequence[str] | None = "f3a9c7d2e481"
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None


def upgrade() -> None:
    """Upgrade schema."""
    op.add_column("submission_states", sa.Column("failure_reason", sa.String(), nullable=True))


def downgrade() -> None:
    """Downgrade schema."""
    raise RuntimeError("Downgrades not supported.")
