"""migrate detailed_qc_results to jsonb

Revision ID: d4e5f6a7b8c9
Revises: fb3df229a77b, 09602efd9105
Create Date: 2026-05-30 14:00:00.000000

"""

from collections.abc import Sequence

import sqlalchemy as sa
from alembic import op
import sqlalchemy.dialects.postgresql as sa_psql

# revision identifiers, used by Alembic.
revision: str = "d4e5f6a7b8c9"
down_revision: str | Sequence[str] | None = ("fb3df229a77b", "09602efd9105")
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None


def upgrade() -> None:
    # 1. Add metrics column as nullable JSONB
    op.add_column(
        "detailed_qc_results",
        sa.Column("metrics", sa.JSON().with_variant(sa_psql.JSONB, "postgresql"), nullable=True)
    )

    # 2. Migrate data
    # We use a raw SQL update to move the data from columns to the JSONB field.
    # PostgreSQL jsonb_build_object is very handy here.
    op.execute(
        """
        UPDATE detailed_qc_results
        SET metrics = jsonb_build_object(
            'percent_bases_above_quality_threshold_minimum_quality', percent_bases_above_quality_threshold_minimum_quality,
            'percent_bases_above_quality_threshold_percent', percent_bases_above_quality_threshold_percent,
            'percent_bases_above_quality_threshold_passed_qc', percent_bases_above_quality_threshold_passed_qc,
            'percent_bases_above_quality_threshold_percent_deviation', percent_bases_above_quality_threshold_percent_deviation,
            'mean_depth_of_coverage', mean_depth_of_coverage,
            'mean_depth_of_coverage_passed_qc', mean_depth_of_coverage_passed_qc,
            'mean_depth_of_coverage_percent_deviation', mean_depth_of_coverage_percent_deviation,
            'targeted_regions_min_coverage', targeted_regions_min_coverage,
            'targeted_regions_above_min_coverage', targeted_regions_above_min_coverage,
            'targeted_regions_above_min_coverage_passed_qc', targeted_regions_above_min_coverage_passed_qc,
            'targeted_regions_above_min_coverage_percent_deviation', targeted_regions_above_min_coverage_percent_deviation
        )
        """
    )

    # 3. Make metrics NOT NULL
    op.alter_column("detailed_qc_results", "metrics", nullable=False)

    # 4. Drop old columns
    with op.batch_alter_table("detailed_qc_results") as batch_op:
        batch_op.drop_column("percent_bases_above_quality_threshold_minimum_quality")
        batch_op.drop_column("percent_bases_above_quality_threshold_percent")
        batch_op.drop_column("percent_bases_above_quality_threshold_passed_qc")
        batch_op.drop_column("percent_bases_above_quality_threshold_percent_deviation")
        batch_op.drop_column("mean_depth_of_coverage")
        batch_op.drop_column("mean_depth_of_coverage_passed_qc")
        batch_op.drop_column("mean_depth_of_coverage_percent_deviation")
        batch_op.drop_column("targeted_regions_min_coverage")
        batch_op.drop_column("targeted_regions_above_min_coverage")
        batch_op.drop_column("targeted_regions_above_min_coverage_passed_qc")
        batch_op.drop_column("targeted_regions_above_min_coverage_percent_deviation")


def downgrade() -> None:
    # 1. Add columns back as nullable
    with op.batch_alter_table("detailed_qc_results") as batch_op:
        batch_op.add_column(sa.Column("percent_bases_above_quality_threshold_minimum_quality", sa.Float(), nullable=True))
        batch_op.add_column(sa.Column("percent_bases_above_quality_threshold_percent", sa.Float(), nullable=True))
        batch_op.add_column(sa.Column("percent_bases_above_quality_threshold_passed_qc", sa.Boolean(), nullable=True))
        batch_op.add_column(sa.Column("percent_bases_above_quality_threshold_percent_deviation", sa.Float(), nullable=True))
        batch_op.add_column(sa.Column("mean_depth_of_coverage", sa.Float(), nullable=True))
        batch_op.add_column(sa.Column("mean_depth_of_coverage_passed_qc", sa.Boolean(), nullable=True))
        batch_op.add_column(sa.Column("mean_depth_of_coverage_percent_deviation", sa.Float(), nullable=True))
        batch_op.add_column(sa.Column("targeted_regions_min_coverage", sa.Float(), nullable=True))
        batch_op.add_column(sa.Column("targeted_regions_above_min_coverage", sa.Float(), nullable=True))
        batch_op.add_column(sa.Column("targeted_regions_above_min_coverage_passed_qc", sa.Boolean(), nullable=True))
        batch_op.add_column(sa.Column("targeted_regions_above_min_coverage_percent_deviation", sa.Float(), nullable=True))

    # 2. Migrate data back
    op.execute(
        """
        UPDATE detailed_qc_results
        SET 
            percent_bases_above_quality_threshold_minimum_quality = (metrics->>'percent_bases_above_quality_threshold_minimum_quality')::float,
            percent_bases_above_quality_threshold_percent = (metrics->>'percent_bases_above_quality_threshold_percent')::float,
            percent_bases_above_quality_threshold_passed_qc = (metrics->>'percent_bases_above_quality_threshold_passed_qc')::boolean,
            percent_bases_above_quality_threshold_percent_deviation = (metrics->>'percent_bases_above_quality_threshold_percent_deviation')::float,
            mean_depth_of_coverage = (metrics->>'mean_depth_of_coverage')::float,
            mean_depth_of_coverage_passed_qc = (metrics->>'mean_depth_of_coverage_passed_qc')::boolean,
            mean_depth_of_coverage_percent_deviation = (metrics->>'mean_depth_of_coverage_percent_deviation')::float,
            targeted_regions_min_coverage = (metrics->>'targeted_regions_min_coverage')::float,
            targeted_regions_above_min_coverage = (metrics->>'targeted_regions_above_min_coverage')::float,
            targeted_regions_above_min_coverage_passed_qc = (metrics->>'targeted_regions_above_min_coverage_passed_qc')::boolean,
            targeted_regions_above_min_coverage_percent_deviation = (metrics->>'targeted_regions_above_min_coverage_percent_deviation')::float
        """
    )

    # 3. Make columns NOT NULL
    with op.batch_alter_table("detailed_qc_results") as batch_op:
        batch_op.alter_column("percent_bases_above_quality_threshold_minimum_quality", nullable=False)
        batch_op.alter_column("percent_bases_above_quality_threshold_percent", nullable=False)
        batch_op.alter_column("percent_bases_above_quality_threshold_passed_qc", nullable=False)
        batch_op.alter_column("percent_bases_above_quality_threshold_percent_deviation", nullable=False)
        batch_op.alter_column("mean_depth_of_coverage", nullable=False)
        batch_op.alter_column("mean_depth_of_coverage_passed_qc", nullable=False)
        batch_op.alter_column("mean_depth_of_coverage_percent_deviation", nullable=False)
        batch_op.alter_column("targeted_regions_min_coverage", nullable=False)
        batch_op.alter_column("targeted_regions_above_min_coverage", nullable=False)
        batch_op.alter_column("targeted_regions_above_min_coverage_passed_qc", nullable=False)
        batch_op.alter_column("targeted_regions_above_min_coverage_percent_deviation", nullable=False)

    # 4. Drop metrics column
    op.drop_column("detailed_qc_results", "metrics")
