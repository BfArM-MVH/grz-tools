"""Command for generating various reports related to GRZ activities."""

import calendar
import csv
import datetime
import itertools
import logging
import re
import typing
from collections import defaultdict
from operator import attrgetter
from pathlib import Path

import click
import sqlalchemy as sa
from grz_common.cli import config_file
from grz_db.models.submission import (
    ChangeRequestEnum,
    ChangeRequestLog,
    ConsentRecord,
    DetailedQCResult,
    Submission,
    SubmissionDb,
    SubmissionStateEnum,
)
from grz_pydantic_models.submission.metadata import GenomicStudyType, Relation, SubmissionType
from sqlalchemy import func as sqlfn
from sqlmodel import select

from ..models.config import DbConfig
from .db.cli import get_submission_db_instance

log = logging.getLogger(__name__)


@click.group()
@config_file
@click.pass_context
def report(ctx: click.Context, config_file: str):
    """
    Generate various reports related to GRZ activities.
    """
    config = DbConfig.from_path(config_file).db
    if not config:
        raise ValueError("DB config not found")

    ctx.obj = {"db_url": config.database_url}


@report.command()
@click.option(
    "--since",
    "since",
    type=datetime.date.fromisoformat,
    help="First date on which to include submissions (default: a week before 'until').",
)
@click.option(
    "--until",
    "until",
    type=datetime.date.fromisoformat,
    help="Last date on which to include submissions (default: today).",
)
@click.option("-s", "--separator", type=str, default="\t", help="Separator between columns (default: tab).")
@click.pass_context
def processed(ctx: click.Context, since: datetime.date | None, until: datetime.date | None, separator: str):
    """
    Generate a report of processed submissions.
    Generally, this is for regular reporting to LEs.
    """
    db = ctx.obj["db_url"]
    submission_db = get_submission_db_instance(db)

    if until is None:
        # default to today
        until = datetime.date.today()

    if since is None:
        # default to a week before 'until'
        since = until - datetime.timedelta(weeks=1)

    submissions = submission_db.list_processed_between(start=since, end=until)

    status_map: dict[bool | None, str] = {
        True: "yes",
        False: "no",
        None: "",
    }

    click.echo(f"# Submissions processed between {since} and {until}")
    click.echo(separator.join(["Submission ID", "Basic QC Passed", "Detailed QC Passed", "Prüfbericht Submitted"]))
    for submission in submissions:
        last_reported_state_change = submission.get_latest_state(filter_to_type=SubmissionStateEnum.REPORTED)
        if not last_reported_state_change:
            continue
        click.echo(
            separator.join(
                [
                    submission.id,
                    status_map[submission.basic_qc_passed],
                    status_map[submission.detailed_qc_passed],
                    last_reported_state_change.timestamp.date().isoformat(),
                ]
            )
        )


def _get_quarter_date_bounds(year: int, quarter: int) -> tuple[datetime.date, datetime.date]:
    quarter_start_date = datetime.date(year=year, month=((quarter - 1) * 3) + 1, day=1)
    quarter_end_month = quarter_start_date.month + 2
    _quarter_end_month_first_weekday, days_in_quarter_end_month = calendar.monthrange(year, quarter_end_month)
    quarter_end_date = quarter_start_date.replace(month=quarter_end_month, day=days_in_quarter_end_month)

    return quarter_start_date, quarter_end_date


def _get_consent_revocations(
    session, quarter_start_date, quarter_end_date
) -> defaultdict[tuple[str, str, str, str], int]:
    """
    Revocation is defined as consent state changing from True to False for a
    donor based on metadata in a submission from this reporting quarter.

    A donor pseudonym is only unique within a single submitter.

    The basic algorithm is:
    1.) Build a map of donors to their consent state at the end of the current quarter.
    2.) Build a map of donors to their consent state at the end of the prior quarter.
    3.) Return donors that were consented at end of prior quarter but are no longer as of end of current quarter.
    """
    subquery_quarter_submissions = (
        select(Submission)
        .where(Submission.submission_date.between(quarter_start_date, quarter_end_date))  # type: ignore[union-attr]
        .subquery()
    )
    query_quarter_consent_records = (
        select(
            subquery_quarter_submissions.c.data_node_id,
            subquery_quarter_submissions.c.submitter_id,
            subquery_quarter_submissions.c.pseudonym,
            ConsentRecord,
        )
        .join(subquery_quarter_submissions, subquery_quarter_submissions.c.id == ConsentRecord.submission_id)
        .order_by(subquery_quarter_submissions.c.submission_date)
    )
    quarter_consent_records = session.exec(query_quarter_consent_records).all()

    # key is (data_node_id, submitter_id, pseudonym, mv|research)
    end_of_quarter_consent_state_by_donor: dict[tuple[str, str, str, str], bool] = {}
    for data_node_id, submitter_id, index_pseudonym, record in quarter_consent_records:
        # iterate over consent records from earliest submission to latest, later overriding earlier
        end_of_quarter_consent_state_by_donor[
            (
                data_node_id,
                submitter_id,
                index_pseudonym if record.relation == Relation.index_ else record.pseudonym,
                "mv",
            )
        ] = record.mv_consented
        end_of_quarter_consent_state_by_donor[
            (
                data_node_id,
                submitter_id,
                index_pseudonym if record.relation == Relation.index_ else record.pseudonym,
                "research",
            )
        ] = record.research_consented

    # now query before the current quarter
    subquery_prior_submissions = select(Submission).where(Submission.submission_date < quarter_start_date).subquery()
    query_prior_consent_records = (
        select(
            subquery_prior_submissions.c.data_node_id,
            subquery_prior_submissions.c.submitter_id,
            subquery_prior_submissions.c.pseudonym,
            ConsentRecord,
        )
        .join(subquery_prior_submissions, subquery_prior_submissions.c.id == ConsentRecord.submission_id)
        .order_by(subquery_prior_submissions.c.submission_date)
    )
    prior_consent_records = session.exec(query_prior_consent_records).all()

    prior_consent_state_by_donor: dict[tuple[str, str, str, str], bool] = {}
    for data_node_id, submitter_id, index_pseudonym, record in prior_consent_records:
        # iterate over consent records from earliest submission to latest, later overriding earlier
        prior_consent_state_by_donor[
            (
                data_node_id,
                submitter_id,
                index_pseudonym if record.relation == Relation.index_ else record.pseudonym,
                "mv",
            )
        ] = record.mv_consented
        prior_consent_state_by_donor[
            (
                data_node_id,
                submitter_id,
                index_pseudonym if record.relation == Relation.index_ else record.pseudonym,
                "research",
            )
        ] = record.research_consented

    number_of_consent_revocations: defaultdict[tuple[str, str, str, str], int] = defaultdict(int)
    for (
        data_node_id,
        submitter_id,
        pseudonym,
        consent_type,
    ), consented in end_of_quarter_consent_state_by_donor.items():
        previously_consented = prior_consent_state_by_donor.get(
            (data_node_id, submitter_id, pseudonym, consent_type), False
        )
        if previously_consented and not consented:
            number_of_consent_revocations[(data_node_id, submitter_id, pseudonym, consent_type)] += 1

    return number_of_consent_revocations


def _dump_overview_report(output_path: Path, database: SubmissionDb, year: int, quarter: int) -> None:
    quarter_start_date, quarter_end_date = _get_quarter_date_bounds(year=year, quarter=quarter)

    node_submitter_id_combos: set[tuple[str | None, str | None]] = set()
    with database._get_session() as session:
        # number_of_end-to-end_tests
        stmt_number_of_end_to_end_tests = (
            select(Submission.data_node_id, Submission.submitter_id, sqlfn.count(1))
            .where(Submission.submission_date.between(quarter_start_date, quarter_end_date))  # type: ignore[union-attr]
            .where(Submission.submission_type == SubmissionType.test)
            .group_by(Submission.data_node_id, Submission.submitter_id)  # type: ignore[arg-type]
        )
        number_of_end_to_end_tests = {
            (node, submitter): count for node, submitter, count in session.exec(stmt_number_of_end_to_end_tests).all()
        }
        node_submitter_id_combos.update(number_of_end_to_end_tests.keys())

        # number_of_passed_end-to-end_tests
        stmt_number_of_passed_end_to_end_tests = (
            select(Submission.data_node_id, Submission.submitter_id, sqlfn.count(1))
            .where(Submission.submission_date.between(quarter_start_date, quarter_end_date))  # type: ignore[union-attr]
            .where(Submission.submission_type == SubmissionType.test)
            .filter(Submission.basic_qc_passed)  # type: ignore[arg-type]
            .group_by(Submission.data_node_id, Submission.submitter_id)  # type: ignore[arg-type]
        )
        number_of_passed_end_to_end_tests = {
            (node, submitter): count
            for node, submitter, count in session.execute(stmt_number_of_passed_end_to_end_tests).all()
        }
        node_submitter_id_combos.update(number_of_passed_end_to_end_tests.keys())

        # number_of_submissions_*
        stmt_number_of_submissions = (
            select(Submission.data_node_id, Submission.submitter_id, Submission.genomic_study_type, sqlfn.count(1))
            .where(Submission.submission_date.between(quarter_start_date, quarter_end_date))  # type: ignore[union-attr]
            .group_by(Submission.data_node_id, Submission.submitter_id, Submission.genomic_study_type)  # type: ignore[arg-type]
        )
        number_of_submissions_by_genomic_study_type = {
            (node, submitter, study_type): count
            for node, submitter, study_type, count in session.execute(stmt_number_of_submissions).all()
        }
        number_of_submissions: defaultdict[tuple[str, str], int] = defaultdict(int)
        for (node, submitter, _), count in number_of_submissions_by_genomic_study_type.items():
            number_of_submissions[(node, submitter)] += count
            node_submitter_id_combos.add((node, submitter))

        # number_of_failed_qcs
        stmt_number_of_failed_qcs = (
            select(Submission.data_node_id, Submission.submitter_id, sqlfn.count(1))
            .where(Submission.submission_date.between(quarter_start_date, quarter_end_date))  # type: ignore[union-attr]
            .filter(sa.not_(Submission.detailed_qc_passed))  # type: ignore[call-overload]
            .group_by(Submission.data_node_id, Submission.submitter_id)  # type: ignore[arg-type]
        )
        number_of_failed_qcs = {
            (node, submitter): count for node, submitter, count in session.execute(stmt_number_of_failed_qcs).all()
        }
        node_submitter_id_combos.update(number_of_failed_qcs.keys())

        # number_of_*_consent_revocations_*
        number_of_consent_revocations = _get_consent_revocations(session, quarter_start_date, quarter_end_date)
        for node, submitter, _, _ in number_of_consent_revocations:
            node_submitter_id_combos.add((node, submitter))

        # number_of_deletions
        stmt_number_of_deletions = (
            select(Submission.data_node_id, Submission.submitter_id, sqlfn.count(1))
            .where(Submission.submission_date.between(quarter_start_date, quarter_end_date))  # type: ignore[union-attr]
            .join(
                select(ChangeRequestLog.submission_id)
                .where(ChangeRequestLog.change == ChangeRequestEnum.DELETE)
                .subquery()
            )
            .group_by(Submission.data_node_id, Submission.submitter_id)  # type: ignore[arg-type]
        )
        number_of_deletions = {
            (node, submitter): count for node, submitter, count in session.execute(stmt_number_of_deletions).all()
        }
        node_submitter_id_combos.update(number_of_deletions.keys())

    with open(output_path, mode="w", encoding="utf-8", newline="") as output_file:
        writer = csv.writer(output_file, delimiter="\t")
        # header
        writer.writerow(
            [
                "genomicDataCenterId",
                "quarter",
                "year",
                "submitterId",
                "number_of_end-to-end_tests",
                "number_of_passed_end-to-end_tests",
                "number_of_submissions_total",
                "number_of_submissions_single",
                "number_of_submissions_duo",
                "number_of_submissions_trio",
                "number_of_failed_qcs",
                "number_of_mv_consent_revocations_index",
                "number_of_research_consent_revocations_index",
                "number_of_mv_consent_revocations_not_index",
                "number_of_research_consent_revocations_not_index",
                "number_of_deletions",
            ]
        )
        for data_node_id, submitter_id in sorted(node_submitter_id_combos):
            if data_node_id is None:
                raise ValueError("At least one submission for the reporting quarter has a null data_node_id")
            if submitter_id is None:
                raise ValueError("At least one submission for the reporting quarter has a null submitter_id")

            writer.writerow(
                [
                    data_node_id,
                    quarter,
                    year,
                    submitter_id,
                    number_of_end_to_end_tests.get((data_node_id, submitter_id), 0),
                    number_of_passed_end_to_end_tests.get((data_node_id, submitter_id), 0),
                    number_of_submissions[(data_node_id, submitter_id)],
                    number_of_submissions_by_genomic_study_type.get(
                        (data_node_id, submitter_id, GenomicStudyType.single), 0
                    ),
                    number_of_submissions_by_genomic_study_type.get(
                        (data_node_id, submitter_id, GenomicStudyType.duo), 0
                    ),
                    number_of_submissions_by_genomic_study_type.get(
                        (data_node_id, submitter_id, GenomicStudyType.trio), 0
                    ),
                    number_of_failed_qcs.get((data_node_id, submitter_id), 0),
                    number_of_consent_revocations.get((data_node_id, submitter_id, "mv", "index"), 0),
                    number_of_consent_revocations.get((data_node_id, submitter_id, "research", "index"), 0),
                    number_of_consent_revocations.get((data_node_id, submitter_id, "mv", "not-index"), 0),
                    number_of_consent_revocations.get((data_node_id, submitter_id, "research", "not-index"), 0),
                    number_of_deletions.get((data_node_id, submitter_id), 0),
                ]
            )


def _dump_dataset_report(output_path: Path, database: SubmissionDb, year: int, quarter: int) -> None:
    quarter_start_date, quarter_end_date = _get_quarter_date_bounds(year=year, quarter=quarter)

    with database._get_session() as session:
        query_quarter_submissions = (
            select(Submission).where(Submission.submission_date.between(quarter_start_date, quarter_end_date))  # type: ignore[union-attr]
        )
        submissions = session.exec(query_quarter_submissions).all()

        subquery_quarter_submissions = query_quarter_submissions.subquery()
        query_consent_records = select(ConsentRecord).join(
            subquery_quarter_submissions, subquery_quarter_submissions.c.id == ConsentRecord.submission_id
        )
        consent_records = session.exec(query_consent_records).all()
        id2no_scope_justifications = {}
        for submission_id, submission_consent_records in itertools.groupby(
            sorted(consent_records, key=attrgetter("submission_id")), key=attrgetter("submission_id")
        ):
            justifications = []
            for record in submission_consent_records:
                justifications.append(
                    "NA"
                    if record.research_consent_missing_justification is None
                    else record.research_consent_missing_justification
                )
            id2no_scope_justifications[submission_id] = ";".join(justifications)

    with open(output_path, mode="w", encoding="utf-8", newline="") as output_file:
        writer = csv.writer(output_file, delimiter="\t")
        # header
        writer.writerow(
            [
                "genomicDataCenterId",
                "quarter",
                "year",
                "submitterId",
                "submissionType",
                "coverageType",
                "diseaseType",
                "sequenceType",
                "libraryType",
                "data_quality_check_passed",
                "has_mv_consent",
                "has_research_consent",
                "researchConsent[].NoScopeJustification",
                "detailed_qc_passed",
                "genomicStudyType",
                "genomicStudySubtype",
                "relation",
                "sequenceSubtype",
            ]
        )
        for submission in submissions:
            # TODO: make this an enum
            detailed_qc_passed = "not_performed"
            if submission.detailed_qc_passed is not None:
                detailed_qc_passed = "yes" if submission.detailed_qc_passed else "no"

            writer.writerow(
                [
                    submission.data_node_id,
                    quarter,
                    year,
                    submission.submitter_id,
                    submission.submission_type,
                    submission.coverage_type,
                    submission.disease_type,
                    ";".join(sorted(submission.sequence_types_index))
                    if submission.sequence_types_index is not None
                    else "NA",
                    ";".join(sorted(submission.library_types_index))
                    if submission.library_types_index is not None
                    else "NA",
                    "yes" if submission.basic_qc_passed else "no",
                    "yes",  # currently MV consent is required to validate, need to update this
                    "yes" if submission.consented else "no",
                    id2no_scope_justifications[submission.id],
                    detailed_qc_passed,
                    submission.genomic_study_type,
                    submission.genomic_study_subtype,
                    "index",
                    ";".join(sorted(submission.sequence_subtypes_index))
                    if submission.sequence_subtypes_index is not None
                    else "NA",
                ]
            )


def _dump_qc_report(output_path: Path, database: SubmissionDb, year: int, quarter: int) -> None:
    quarter_start_date, quarter_end_date = _get_quarter_date_bounds(year=year, quarter=quarter)

    with database._get_session() as session:
        query_submissions_that_failed_detailed_qc = (
            select(Submission)
            .where(Submission.submission_date.between(quarter_start_date, quarter_end_date))  # type: ignore[union-attr]
            .filter(sa.not_(Submission.detailed_qc_passed))  # type: ignore[call-overload]
        )
        submissions_that_failed_detailed_qc = session.exec(query_submissions_that_failed_detailed_qc).all()
        query_reports_of_failed_submissions = select(DetailedQCResult).join(
            query_submissions_that_failed_detailed_qc.subquery()
        )
        reports_of_failed_submissions = session.exec(query_reports_of_failed_submissions).all()

    id2submission = {submission.id: submission for submission in submissions_that_failed_detailed_qc}
    with open(output_path, mode="w", encoding="utf-8", newline="") as output_file:
        writer = csv.writer(output_file, delimiter="\t")
        # header
        writer.writerow(
            [
                "genomicDataCenterId",
                "quarter",
                "year",
                "submitterId",
                "submissionDate",
                "submissionType",
                "diseaseType",
                "genomicStudyType",
                "genomicStudySubtype",
                "relation",
                "sequenceType",
                "sequenceSubtype",
                "libraryType",
                "percentBasesAboveQualityThreshold.minimumQuality",
                "percentBasesAboveQualityThreshold.percent",
                "percentBasesAboveQualityThreshold_detailedQC_passed",
                "percentBasesAboveQualityThreshold_detailedQC_deviation%",
                "meanDepthOfCoverage",
                "meanDepthOfCoverage_detailedQC_passed",
                "meanDepthOfCoverage_detailedQC_deviation%",
                "minCoverage",
                "targetedRegionsAboveMinCoverage",
                "targetedRegionsAboveMinCoverage_detailedQC_passed",
                "targetedRegionsAboveMinCoverage_detailedQC_deviation%",
            ]
        )

        for report in reports_of_failed_submissions:
            submission = id2submission[report.submission_id]
            # TODO: join on consent table with donor pseudonym to get relation
            lab_datum_id_match = re.match(r"([A-z]+)([0-9]+)_([A-z]+)([0-9]+)", report.lab_datum_id)
            if lab_datum_id_match is None:
                raise ValueError(f"Failed to parse lab datum ID '{report.lab_datum_id}'")
            relation = Relation(lab_datum_id_match.group(1))
            writer.writerow(
                [
                    submission.data_node_id,
                    quarter,
                    year,
                    submission.submitter_id,
                    submission.submission_date,
                    submission.submission_type,
                    submission.disease_type,
                    submission.genomic_study_type,
                    submission.genomic_study_subtype,
                    relation,
                    report.sequence_type,
                    report.sequence_subtype,
                    report.library_type,
                    report.percent_bases_above_quality_threshold_minimum_quality,
                    report.percent_bases_above_quality_threshold_percent,
                    "yes" if report.percent_bases_above_quality_threshold_passed_qc else "no",
                    report.percent_bases_above_quality_threshold_percent_deviation,
                    report.mean_depth_of_coverage,
                    "yes" if report.mean_depth_of_coverage_passed_qc else "no",
                    report.mean_depth_of_coverage_percent_deviation,
                    report.targeted_regions_min_coverage,
                    report.targeted_regions_above_min_coverage,
                    "yes" if report.targeted_regions_above_min_coverage_passed_qc else "no",
                    report.targeted_regions_above_min_coverage_percent_deviation,
                ]
            )


@report.command()
@click.option(
    "--quarter",
    "quarter",
    type=click.IntRange(min=1, max=4),
    help="Quarter to generate report for.",
)
@click.option(
    "--year",
    "year",
    type=click.IntRange(min=2025),
    help="Year to generate report for.",
)
@click.option(
    "--outdir",
    "output_directory",
    type=click.Path(
        exists=True, file_okay=False, dir_okay=True, readable=True, writable=True, resolve_path=True, path_type=Path
    ),
    default=Path.cwd(),
    help="Directory to output TSV files. Defaults to current directory.",
)
@click.pass_context
def quarterly(ctx: click.Context, year: int | None, quarter: int | None, output_directory: Path):
    """
    Generate the tables for the quarterly report.
    """
    db = ctx.obj["db_url"]
    submission_db = get_submission_db_instance(db)

    if bool(year) != bool(quarter):
        raise click.UsageError("Both year and quarter must be provided or omitted.")

    if (year and quarter) is None:
        today = datetime.date.today()
        quarter = ((today.month - 1) % 3) + 1
        # default to last quarter if ended less than 15 days ago otherwise current quarter
        if today <= datetime.date(year=today.year, month=1, day=15):
            year = today.year - 1
            quarter = 4
        else:
            year = today.year

        if (today.month in {4, 7, 10}) and (today.day <= 15):
            quarter -= 1

    # help out the type checker
    year = typing.cast(int, year)
    quarter = typing.cast(int, quarter)

    log.info("Generating quarterly report for Q%d %d", quarter, year)

    overview_output_path = output_directory / "Gesamtübersicht.tsv"
    _dump_overview_report(overview_output_path, submission_db, year, quarter)

    dataset_output_path = output_directory / "Infos_zu_Datensätzen.tsv"
    _dump_dataset_report(dataset_output_path, submission_db, year, quarter)

    qc_output_path = output_directory / "Detailprüfung.tsv"
    _dump_qc_report(qc_output_path, submission_db, year, quarter)
