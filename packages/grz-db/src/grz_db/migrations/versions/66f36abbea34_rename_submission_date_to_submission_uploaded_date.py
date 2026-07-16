"""rename submission_date to submission_uploaded_date

Revision ID: 66f36abbea34
Revises: 09602efd9105
Create Date: 2026-06-10 18:20:26.144612+00:00

"""

from collections.abc import Sequence

from alembic import op

# revision identifiers, used by Alembic.
revision: str = "66f36abbea34"
down_revision: str | Sequence[str] | None = "09602efd9105"
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None


def upgrade() -> None:
    """Upgrade schema."""
    op.alter_column("submissions", "submission_date", new_column_name="submission_uploaded_date")


def downgrade() -> None:
    """Downgrade schema."""
    raise RuntimeError("Downgrades not supported.")
