"""add grzctl_version and qc_workflow_version columns

Revision ID: c5f6a8b9d2e1
Revises: 834bd50b8734
Create Date: 2026-05-06 00:00:00.000000
"""

from collections.abc import Sequence

import sqlalchemy as sa
from alembic import op

revision: str = "c5f6a8b9d2e1"
down_revision: str | Sequence[str] | None = "834bd50b8734"
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None


def upgrade() -> None:
    # Add grzctl_versions to submission_states
    op.add_column(
        "submission_states",
        sa.Column("grzctl_versions", sa.JSON(), nullable=True),
    )

    # Add qc_workflow_version to detailed_qc_results
    op.add_column(
        "detailed_qc_results",
        sa.Column("qc_workflow_version", sa.String(length=64), nullable=True),
    )


def downgrade() -> None:
    """Downgrade schema."""
    raise RuntimeError("Downgrades not supported.")
