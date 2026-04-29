"""merge_main_and_failure_reason_branch

Revision ID: 09602efd9105
Revises: 834bd50b8734, a4de712d8faa
Create Date: 2026-04-29 09:35:54.981907+00:00

"""

from collections.abc import Sequence

# revision identifiers, used by Alembic.
revision: str = "09602efd9105"
down_revision: str | Sequence[str] | None = ("834bd50b8734", "a4de712d8faa")
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None


def upgrade() -> None:
    """Upgrade schema."""
    pass


def downgrade() -> None:
    """Downgrade schema."""
    raise RuntimeError("Downgrades not supported.")
