"""add qc_workflow_version to detailed_qc_results

Revision ID: b4e8a1c9d2f7
Revises: a9f3c7e2b6d1
Create Date: 2026-05-06 00:00:00.000000
"""

from collections.abc import Sequence

import sqlalchemy as sa
from alembic import op

revision: str = "b4e8a1c9d2f7"
down_revision: str | Sequence[str] | None = "a9f3c7e2b6d1"
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None


def upgrade() -> None:
    op.add_column(
        "detailed_qc_results",
        sa.Column("qc_workflow_version", sa.String(length=64), nullable=True),
    )


def downgrade() -> None:
    op.drop_column("detailed_qc_results", "qc_workflow_version")
