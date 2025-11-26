"""initial Quartalsbericht support

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
    # new donors table
    op.create_table(
        "donors",
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
            primary_key=True,
            nullable=False,
        ),
        sa.Column(
            "relation",
            sa.Enum("mother", "father", "brother", "sister", "child", "index_", "other", name="relation"),
            nullable=False,
        ),
        sa.Column("library_types", AutoString(), nullable=False),
        sa.Column("sequence_types", AutoString(), nullable=False),
        sa.Column("sequence_subtypes", AutoString(), nullable=False),
        sa.Column("mv_consented", sa.Boolean(), nullable=False),
        sa.Column("research_consented", sa.Boolean(), nullable=False),
        sa.Column("research_consent_missing_justifications", AutoString(), nullable=True),
    )

    # explicitly create enums if needed (PostgreSQL)
    genomic_study_type_enum = sa.Enum("single", "duo", "trio", name="genomicstudytype")
    genomic_study_type_enum.create(op.get_bind())
    genomic_study_subtype_enum = sa.Enum("tumor_only", "tumor_germline", "germline_only", name="genomicstudysubtype")
    genomic_study_subtype_enum.create(op.get_bind())
    coverage_type_enum = sa.Enum(
        "GKV", "PKV", "BG", "SEL", "SOZ", "GPV", "PPV", "BEI", "SKT", "UNK", name="coveragetype"
    )
    coverage_type_enum.create(op.get_bind())

    # we have to batch alter because SQLite doesn't support ALTER on column types
    with op.batch_alter_table("submissions") as batch_op:
        # will need to be repopulated in donors table before reporting, so just drop
        batch_op.drop_column("library_type")
        # modified/extra submission table columns
        batch_op.add_column(sa.Column("genomic_study_type", genomic_study_type_enum))
        batch_op.add_column(
            sa.Column(
                "genomic_study_subtype",
                genomic_study_subtype_enum,
            ),
        )
        batch_op.add_column(sa.Column("coverage_type", coverage_type_enum))

    # already should exist so don't need to re-create, just define
    # but must be postgresql.ENUM until fixed in sqlalchemy 2.1
    # see https://github.com/sqlalchemy/alembic/issues/1347
    library_type_enum = sa.dialects.postgresql.ENUM(
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
        create_type=False,
    )

    # new detailed QC results table
    op.create_table(
        "detailed_qc_results",
        sa.Column(
            "submission_id",
            AutoString(),
            primary_key=True,
            nullable=False,
        ),
        sa.Column("lab_datum_id", AutoString(), primary_key=True, nullable=False),
        sa.Column("pseudonym", AutoString(), nullable=False),
        sa.Column("timestamp", sa.DateTime(timezone=True), nullable=False, primary_key=True),
        sa.Column("sequence_type", sa.Enum("dna", "rna", name="sequencetype"), nullable=False),
        sa.Column(
            "sequence_subtype",
            sa.Enum("germline", "somatic", "other", "unknown", name="sequencesubtype"),
            nullable=False,
        ),
        sa.Column(
            "library_type",
            library_type_enum,
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
        sa.ForeignKeyConstraint(["submission_id", "pseudonym"], ["donors.submission_id", "donors.pseudonym"]),
    )


def downgrade() -> None:
    """Downgrade schema."""
    raise RuntimeError("Downgrades not supported.")
