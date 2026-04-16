"""Add failure_reason column to submission_states

Revision ID: c063b408d262
Revises: 8aa6fc8d118a
Create Date: 2026-04-13 09:51:49.819155+00:00

"""

from collections.abc import Sequence

import sqlalchemy as sa
from alembic import op
from grz_db.models.submission import FailureReasonEnum

# revision identifiers, used by Alembic.
revision: str = "c063b408d262"
down_revision: str | Sequence[str] | None = "8aa6fc8d118a"
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None


def upgrade() -> None:
    """Upgrade schema - added the failure_reason column to submission_states."""
    failure_reason_enum = sa.Enum(*[e.value for e in FailureReasonEnum], name="failure_reason_enum")
    failure_reason_enum.create(op.get_bind())
    op.add_column("submission_states", sa.Column("failure_reason", failure_reason_enum, nullable=True))


def downgrade() -> None:
    """Downgrade schema."""
    raise RuntimeError("Downgrades not supported.")
