"""add_reporting_columns

Revision ID: 9bf36a91c87e
Revises:
Create Date: 2025-07-24 09:28:26.368617+00:00

"""

from collections.abc import Sequence

import sqlalchemy as sa
from alembic import op
from sqlmodel.sql.sqltypes import AutoString

# revision identifiers, used by Alembic.
revision: str = "9bf36a91c87e"
down_revision: str | Sequence[str] | None = "1a9bd994df1b"
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None


def upgrade() -> None:
    """Upgrade schema."""
    op.add_column("submissions", sa.Column("submission_date", sa.Date(), nullable=True))
    op.add_column(
        "submissions",
        sa.Column(
            "submission_type", sa.Enum("initial", "followup", "addition", "correction", "test", name="submissiontype")
        ),
    )
    op.add_column("submissions", sa.Column("submitter_id", AutoString()))
    op.add_column("submissions", sa.Column("data_node_id", AutoString()))
    op.add_column(
        "submissions", sa.Column("disease_type", sa.Enum("oncological", "rare", "hereditary", name="diseasetype"))
    )
    op.add_column(
        "submissions",
        sa.Column(
            "library_type",
            sa.Enum(
                "panel",
                "panel_lr",
                "wes",
                "wes_lr",
                "wgs",
                "wgs_lr",
                "wxs",
                "wxs_lr",
                "other",
                "unknown",
                name="librarytype",
            ),
        ),
    )
    op.add_column("submissions", sa.Column("basic_qc_passed", sa.Boolean()))
    op.add_column("submissions", sa.Column("consented", sa.Boolean()))
    op.add_column("submissions", sa.Column("detailed_qc_passed", sa.Boolean()))


def downgrade() -> None:
    """Downgrade schema."""
    raise RuntimeError("Downgrades not supported.")
