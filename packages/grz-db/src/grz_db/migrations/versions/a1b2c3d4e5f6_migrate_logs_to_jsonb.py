"""migrate logs to jsonb

Revision ID: a1b2c3d4e5f6
Revises: e1f2a3b4c5d6
Create Date: 2026-06-03 12:00:00.000000

"""

from collections.abc import Sequence

import sqlalchemy as sa
from alembic import op

# revision identifiers, used by Alembic.
revision: str = "a1b2c3d4e5f6"
down_revision: str | Sequence[str] | None = "e1f2a3b4c5d6"
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None


def upgrade() -> None:
    # Migrate submission_states table
    # PostgreSQL specific ALTER TABLE ... USING ...
    op.execute("ALTER TABLE submission_states ALTER COLUMN data TYPE JSONB USING data::jsonb")
    op.execute("ALTER TABLE submission_states ALTER COLUMN grzctl_versions TYPE JSONB USING grzctl_versions::jsonb")
    
    # Migrate submission_change_requests table
    op.execute("ALTER TABLE submission_change_requests ALTER COLUMN data TYPE JSONB USING data::jsonb")


def downgrade() -> None:
    # Revert submission_states table
    op.execute("ALTER TABLE submission_states ALTER COLUMN data TYPE JSON USING data::json")
    op.execute("ALTER TABLE submission_states ALTER COLUMN grzctl_versions TYPE JSON USING grzctl_versions::json")
    
    # Revert submission_change_requests table
    op.execute("ALTER TABLE submission_change_requests ALTER COLUMN data TYPE JSON USING data::json")
