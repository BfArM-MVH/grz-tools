"""inital Quartalsbericht support

Revision ID: fb3df229a77b
Revises: 9bf36a91c87e
Create Date: 2025-09-12 12:47:43.113455+00:00

"""

from collections.abc import Sequence

import sqlalchemy as sa
from alembic import op
from sqlalchemy.schema import ForeignKey
from sqlmodel.sql.sqltypes import AutoString

# revision identifiers, used by Alembic.
revision: str = "fb3df229a77b"
down_revision: str | Sequence[str] | None = "9bf36a91c87e"
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None


def upgrade() -> None:
    """Upgrade schema."""
    # we have to batch alter because SQLite doesn't support ALTER on column types
    with op.batch_alter_table("submissions") as batch_op:
        # modified/extra submission table columns
        batch_op.alter_column("library_type", type_=AutoString(), new_column_name="library_types_index")
        batch_op.add_column(sa.Column("genomic_study_type", sa.Enum("single", "duo", "trio", name="genomicstudytype")))
        batch_op.add_column(
            sa.Column(
                "genomic_study_subtype",
                sa.Enum("tumor_only", "tumor_germline", "germline_only", name="genomicstudysubtype"),
            ),
        )
        batch_op.add_column(
            sa.Column(
                "coverage_type",
                sa.Enum("GKV", "PKV", "BG", "SEL", "SOZ", "GPV", "PPV", "BEI", "SKT", "UNK", name="coveragetype"),
            )
        )
        batch_op.add_column(sa.Column("sequence_types_index", AutoString()))
        batch_op.add_column(sa.Column("sequence_subtypes_index", AutoString()))

    # new consent record table
    op.create_table(
        "consent_records",
        sa.Column(
            "submission_id",
            AutoString(),
            ForeignKey("submissions.id"),
            primary_key=True,
            nullable=False,
        ),
        sa.Column(
            "pseudonym",
            AutoString(),
            ForeignKey("submissions.pseudonym"),
            primary_key=True,
            nullable=False,
        ),
        sa.Column(
            "relation",
            sa.Enum("mother", "father", "brother", "sister", "child", "index_", "other", name="relation"),
            nullable=False,
        ),
        sa.Column("mv_consented", sa.Boolean(), nullable=False),
        sa.Column("research_consented", sa.Boolean(), nullable=False),
        sa.Column(
            "research_consent_missing_justification",
            sa.Enum(
                "UNABLE",
                "REFUSED",
                "NO_RETURN",
                "OTHER",
                "LE_TECH",
                "LE_ORG",
                name="researchconsentnoscopejustification",
            ),
        ),
    )

    # new detailed QC results table
    op.create_table(
        "detailed_qc_results",
        sa.Column(
            "submission_id",
            AutoString(),
            ForeignKey("submissions.id"),
            primary_key=True,
            nullable=False,
        ),
        sa.Column("lab_datum_id", AutoString(), primary_key=True, nullable=False),
        sa.Column("timestamp", sa.DateTime(timezone=True), nullable=False),
        sa.Column("sequence_type", sa.Enum("dna", "rna", name="sequencetype"), nullable=False),
        sa.Column(
            "sequence_subtype",
            sa.Enum("germline", "somatic", "other", "unknown", name="sequencesubtype"),
            nullable=False,
        ),
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
            nullable=False,
        ),
        sa.Column("percent_bases_above_quality_threshold_minimum_quality", sa.Float(), nullable=False),
        sa.Column("percent_bases_above_quality_threshold_percent", sa.Float(), nullable=False),
        sa.Column("percent_bases_above_quality_threshold_passed_qc", sa.Boolean(), nullable=False),
        sa.Column("percent_bases_above_quality_threshold_percent_deviation", sa.Float(), nullable=False),
        sa.Column("mean_depth_of_coverage", sa.Float(), nullable=False),
        sa.Column("mean_depth_of_coverage_passed_qc", sa.Boolean(), nullable=False),
        sa.Column("mean_depth_of_coverage_percent_deviation", sa.Float(), nullable=False),
        sa.Column("targeted_regions_min_coverage", sa.Float(), nullable=False),
        sa.Column("targeted_regions_above_min_coverage", sa.Float(), nullable=False),
        sa.Column("targeted_regions_above_min_coverage_passed_qc", sa.Boolean(), nullable=False),
        sa.Column("targeted_regions_above_min_coverage_percent_deviation", sa.Float(), nullable=False),
    )


def downgrade() -> None:
    """Downgrade schema."""
    raise RuntimeError("Downgrades not supported.")
