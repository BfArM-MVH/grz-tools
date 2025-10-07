"""Command for generating various reports related to GRZ activities."""

import calendar
import csv
import datetime
import itertools
import logging
import typing
from collections import defaultdict
from enum import StrEnum
from operator import attrgetter
from pathlib import Path

import click
import sqlalchemy as sa
from grz_common.cli import config_file
from grz_db.models.submission import (
    ChangeRequestEnum,
    ChangeRequestLog,
    DetailedQCResult,
    Donor,
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

    The basic logic for determining revocation is:
    "unconsented at end of quarter" AND ("consented at end of prior quarter" OR "consented sometime within current quarter")
    """
    subquery_quarter_submissions = (
        select(Submission)
        .where(Submission.submission_date.between(quarter_start_date, quarter_end_date))  # type: ignore[union-attr]
        .subquery()
    )
    query_quarter_donors = (
        select(
            subquery_quarter_submissions.c.data_node_id,
            subquery_quarter_submissions.c.submitter_id,
            subquery_quarter_submissions.c.pseudonym,
            Donor,
        )
        .join(subquery_quarter_submissions, subquery_quarter_submissions.c.id == Donor.submission_id)
        .order_by(subquery_quarter_submissions.c.submission_date)
    )
    quarter_donors = session.exec(query_quarter_donors).all()

    # key is (data_node_id, submitter_id, pseudonym, mv|research)
    quarter_consent_states_by_donor: defaultdict[tuple[str, str, str, str], list[bool]] = defaultdict(list)
    # key is (data_node_id, submiter_id, pseudonym)
    relation_by_donor: dict[tuple[str, str, str], Relation] = {}
    # iterate over donors from earliest submission to latest, tracking all consent states
    for data_node_id, submitter_id, index_pseudonym, donor in quarter_donors:
        pseudonym = index_pseudonym if donor.relation == Relation.index_ else donor.pseudonym
        quarter_consent_states_by_donor[
            (
                data_node_id,
                submitter_id,
                pseudonym,
                "mv",
            )
        ].append(donor.mv_consented)
        quarter_consent_states_by_donor[
            (
                data_node_id,
                submitter_id,
                pseudonym,
                "research",
            )
        ].append(donor.research_consented)
        relation_by_donor[(data_node_id, submitter_id, pseudonym)] = donor.relation

    # now query before the current quarter
    subquery_prior_submissions = select(Submission).where(Submission.submission_date < quarter_start_date).subquery()
    query_prior_donors = (
        select(
            subquery_prior_submissions.c.data_node_id,
            subquery_prior_submissions.c.submitter_id,
            subquery_prior_submissions.c.pseudonym,
            Donor,
        )
        .join(subquery_prior_submissions, subquery_prior_submissions.c.id == Donor.submission_id)
        .order_by(subquery_prior_submissions.c.submission_date)
    )
    prior_donors = session.exec(query_prior_donors).all()

    prior_consent_state_by_donor: dict[tuple[str, str, str, str], bool] = {}
    # iterate over donors from earliest submission to latest, later overriding earlier
    for data_node_id, submitter_id, index_pseudonym, donor in prior_donors:
        pseudonym = index_pseudonym if donor.relation == Relation.index_ else donor.pseudonym
        # skip over donors that don't show up in current quarter
        if (data_node_id, submitter_id, pseudonym, "mv") in quarter_consent_states_by_donor:
            prior_consent_state_by_donor[
                (
                    data_node_id,
                    submitter_id,
                    pseudonym,
                    "mv",
                )
            ] = donor.mv_consented
        if (data_node_id, submitter_id, pseudonym, "research") in quarter_consent_states_by_donor:
            prior_consent_state_by_donor[
                (
                    data_node_id,
                    submitter_id,
                    pseudonym,
                    "research",
                )
            ] = donor.research_consented

    number_of_consent_revocations: defaultdict[tuple[str, str, str, str], int] = defaultdict(int)
    for (
        data_node_id,
        submitter_id,
        pseudonym,
        consent_type,
    ), consent_states in quarter_consent_states_by_donor.items():
        previously_consented = prior_consent_state_by_donor.get(
            (data_node_id, submitter_id, pseudonym, consent_type), False
        )
        consented_at_end_of_quarter = consent_states[-1]
        consented_within_quarter = True in consent_states[:-1]
        if not consented_at_end_of_quarter and (previously_consented or consented_within_quarter):
            relation = relation_by_donor[(data_node_id, submitter_id, pseudonym)]
            number_of_consent_revocations[
                (data_node_id, submitter_id, consent_type, "index" if relation == Relation.index_ else "not-index")
            ] += 1

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


class DetailedQCPassedReportState(StrEnum):
    NOT_PERFORMED = "not_performed"
    YES = "yes"
    NO = "no"


def _dump_dataset_report(output_path: Path, database: SubmissionDb, year: int, quarter: int) -> None:
    quarter_start_date, quarter_end_date = _get_quarter_date_bounds(year=year, quarter=quarter)

    with database._get_session() as session:
        query_quarter_submissions = (
            select(Submission).where(Submission.submission_date.between(quarter_start_date, quarter_end_date))  # type: ignore[union-attr]
        )
        submissions = session.exec(query_quarter_submissions).all()

        subquery_quarter_submissions = query_quarter_submissions.subquery()
        query_donors = select(Donor).join(
            subquery_quarter_submissions, subquery_quarter_submissions.c.id == Donor.submission_id
        )
        donors = session.exec(query_donors).all()

        id2sequence_types_index = {}
        id2sequence_subtypes_index = {}
        id2library_types_index = {}
        id2no_scope_justifications = {}
        id2mv_consented = {}
        id2research_consented = {}
        for submission_id, submission_donors in itertools.groupby(
            sorted(donors, key=attrgetter("submission_id")), key=attrgetter("submission_id")
        ):
            sequence_types_index = "NA"
            sequence_subtypes_index = "NA"
            library_types_index = "NA"
            for donor in submission_donors:
                if donor.relation == Relation.index_:
                    sequence_types_index = ";".join(sorted(donor.sequence_types))
                    sequence_subtypes_index = ";".join(sorted(donor.sequence_subtypes))
                    library_types_index = ";".join(sorted(donor.library_types))
                    if donor.research_consent_missing_justifications:
                        id2no_scope_justifications[submission_id] = ";".join(
                            sorted(donor.research_consent_missing_justifications)
                        )
                    id2mv_consented[submission_id] = donor.mv_consented
                    id2research_consented[submission_id] = donor.research_consented

            id2sequence_types_index[submission_id] = sequence_types_index
            id2sequence_subtypes_index[submission_id] = sequence_subtypes_index
            id2library_types_index[submission_id] = library_types_index

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
            detailed_qc_passed = DetailedQCPassedReportState.NOT_PERFORMED
            if submission.detailed_qc_passed is not None:
                detailed_qc_passed = (
                    DetailedQCPassedReportState.YES if submission.detailed_qc_passed else DetailedQCPassedReportState.NO
                )
            if submission.id in id2mv_consented:
                mv_consented = "yes" if id2mv_consented[submission.id] else "no"
            else:
                mv_consented = "NA"

            if submission.id in id2research_consented:
                research_consented = "yes" if id2research_consented[submission.id] else "no"
            else:
                research_consented = "NA"

            writer.writerow(
                [
                    submission.data_node_id,
                    quarter,
                    year,
                    submission.submitter_id,
                    submission.submission_type,
                    submission.coverage_type,
                    submission.disease_type,
                    id2sequence_types_index.get(submission.id, "NA"),
                    id2library_types_index.get(submission.id, "NA"),
                    "yes" if submission.basic_qc_passed else "no",
                    mv_consented,
                    research_consented,
                    id2no_scope_justifications.get(submission.id, "NA"),
                    detailed_qc_passed,
                    submission.genomic_study_type,
                    submission.genomic_study_subtype,
                    "index",
                    id2sequence_subtypes_index.get(submission.id, "NA"),
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

        subquery_submissions_that_failed_detailed_qc = query_submissions_that_failed_detailed_qc.subquery()
        query_reports_of_failed_submissions = select(DetailedQCResult).join(
            subquery_submissions_that_failed_detailed_qc,
            subquery_submissions_that_failed_detailed_qc.c.id == DetailedQCResult.submission_id,
        )
        reports_of_failed_submissions = session.exec(query_reports_of_failed_submissions).all()

        subquery_reports_of_failed_submissions = query_reports_of_failed_submissions.subquery()
        query_relevant_donors = select(Donor).join(
            subquery_reports_of_failed_submissions,
            (subquery_reports_of_failed_submissions.c.submission_id == Donor.submission_id)  # type: ignore[arg-type]
            & (subquery_reports_of_failed_submissions.c.pseudonym == Donor.pseudonym),
        )
        relevant_donors = session.exec(query_relevant_donors).all()

    # Note: for future, this join only works because *both* DetailedQCResult
    # and Donor store index pseudonym as "index", since the submission table
    # has the index pseudonym
    idpseudo2relation = {(donor.submission_id, donor.pseudonym): donor.relation for donor in relevant_donors}
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
            relation = Relation(idpseudo2relation[(report.submission_id, report.pseudonym)])
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


def date_to_quarter_year(date: datetime.date) -> tuple[int, int]:
    """Return (1-based quarter, year) given a date."""
    return ((date.month - 1) % 3) + 1, date.year


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
    help="Directory to output TSV files. Defaults to current directory.",
)
@click.pass_context
def quarterly(ctx: click.Context, year: int | None, quarter: int | None, output_directory: Path | None):
    """
    Generate the tables for the quarterly report.
    """
    db = ctx.obj["db_url"]
    submission_db = get_submission_db_instance(db)

    if bool(year) != bool(quarter):
        raise click.UsageError("Both year and quarter must be provided or omitted.")

    if (year and quarter) is None:
        # default to last quarter if ended less than 15 days ago otherwise current quarter
        check_date = datetime.date.today() - datetime.timedelta(days=15)
        quarter, year = date_to_quarter_year(check_date)

    # help out the type checker
    year = typing.cast(int, year)
    quarter = typing.cast(int, quarter)

    log.info("Generating quarterly report for Q%d %d", quarter, year)

    if output_directory is None:
        output_directory = Path.cwd()

    overview_output_path = output_directory / "Gesamtübersicht.tsv"
    _dump_overview_report(overview_output_path, submission_db, year, quarter)

    dataset_output_path = output_directory / "Infos_zu_Datensätzen.tsv"
    _dump_dataset_report(dataset_output_path, submission_db, year, quarter)

    qc_output_path = output_directory / "Detailprüfung.tsv"
    _dump_qc_report(qc_output_path, submission_db, year, quarter)
