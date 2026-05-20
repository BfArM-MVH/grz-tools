"""add selected_for_qc column

Revision ID: 9023ab55c239
Revises: 8aa6fc8d118a
Create Date: 2026-03-17 08:33:43.113455+00:00

"""

from collections.abc import Sequence

import sqlalchemy as sa
from alembic import op

revision: str = "9023ab55c239"
down_revision: str | Sequence[str] | None = "8aa6fc8d118a"
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None


def upgrade() -> None:
    """Upgrade schema."""
    op.add_column("submissions", sa.Column("selected_for_qc", sa.Boolean(), nullable=True))
    op.execute(
        sa.text(
            """
            UPDATE submissions
            SET selected_for_qc = TRUE
            WHERE id IN (
                SELECT DISTINCT submission_id
                FROM submission_states
                WHERE state IN ('QCING', 'QCED')
            )
            """
        )
    )


def downgrade() -> None:
    """Downgrade schema."""
    raise RuntimeError("Downgrades not supported.")
