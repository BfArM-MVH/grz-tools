"""Command for generating various reports related to GRZ activities."""

import csv
import datetime
import itertools
import json
import logging
import re
import typing
from collections import defaultdict
from enum import StrEnum
from operator import attrgetter
from pathlib import Path

import click
import requests
import sqlalchemy as sa
from grz_db.models.author import Author
from grz_db.models.submission import (
    ChangeRequestEnum,
    ChangeRequestLog,
    DetailedQCResult,
    Donor,
    Submission,
    SubmissionDb,
    SubmissionStateEnum,
)
from grz_pydantic_models.dates import date_to_quarter_year, quarter_date_bounds
from grz_pydantic_models.submission.metadata import GenomicStudyType, Relation, SubmissionType
from grz_pydantic_models.submission.thresholds import PCT_DEV_CUTOFF
from sqlalchemy import func as sqlfn
from sqlalchemy.orm import selectinload
from sqlmodel import select

from ..commands import grzctl_configuration
from ..models.config import GrzctlConfig

log = logging.getLogger(__name__)


def get_submission_db_instance(db_url: str, author: Author | None = None) -> SubmissionDb:
    """Creates and returns an instance of SubmissionDb."""
    return SubmissionDb(db_url=db_url, author=author)


@click.group()
@grzctl_configuration
@click.pass_context
def report(ctx: click.Context, configuration: GrzctlConfig):
    """
    Generate various reports related to GRZ activities.
    """
    db = configuration.db

    ctx.obj = {
        "db_url": db.database_url,
        "grz_id": configuration.identifiers.grz,
        "le_ids": list(configuration.leistungserbringer),
    }


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


def _get_consent_revocations(
    session, quarter_start_date, quarter_end_date
) -> defaultdict[tuple[str, str, str, str], int]:
    """
    Revocation is defined as consent state changing from True to False for a
    donor based on metadata in a submission from this reporting quarter.

    A donor pseudonym is only unique within a single submitter.

    The index donor has no pseudonym of its own: its donorPseudonym is the literal "index", so
    ``submissions.local_case_id`` stands in for it here. ``Submission.pseudonym``, the case's
    psn, takes over once the local case ID is removed; keying on it while no psn is assigned
    would merge every index donor of a submitter into one. Until then the two index counts are
    only as good as that localCaseId. A localCaseId reused across patients merges those
    patients into a single donor.
    That donor's consent history is then several people's, read in upload order. The rule below
    decides a revocation from that mixture. The cases migration logs these keys once, when
    ``grzctl db upgrade`` runs; case tracking refuses to group their submissions for the same
    reason.

    The basic logic for determining revocation is:
    "unconsented at end of quarter" AND ("consented at end of prior quarter" OR "consented sometime within current quarter")
    """
    query_quarter_donors = (
        select(Submission.data_node_id, Submission.submitter_id, Submission.local_case_id, Donor)
        .join_from(Donor, Submission, Submission.id == Donor.submission_id)  # type: ignore[arg-type]
        .where(Submission.submission_uploaded_date.between(quarter_start_date, quarter_end_date))  # type: ignore[union-attr]
        .order_by(Submission.submission_uploaded_date)  # type: ignore[arg-type]
    )
    quarter_donors = session.exec(query_quarter_donors).all()

    # key is (data_node_id, submitter_id, pseudonym, mv|research)
    quarter_consent_states_by_donor: defaultdict[tuple[str, str, str, str], list[bool]] = defaultdict(list)
    # key is (data_node_id, submiter_id, pseudonym)
    relation_by_donor: dict[tuple[str, str, str], Relation] = {}
    # iterate over donors from earliest submission to latest, tracking all consent states
    for data_node_id, submitter_id, index_local_case_id, donor in quarter_donors:
        pseudonym = index_local_case_id if donor.relation == Relation.index_ else donor.pseudonym
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
    query_prior_donors = (
        select(Submission.data_node_id, Submission.submitter_id, Submission.local_case_id, Donor)
        .join_from(Donor, Submission, Submission.id == Donor.submission_id)  # type: ignore[arg-type]
        .where(Submission.submission_uploaded_date < quarter_start_date)
        .order_by(Submission.submission_uploaded_date)  # type: ignore[arg-type]
    )
    prior_donors = session.exec(query_prior_donors).all()

    prior_consent_state_by_donor: dict[tuple[str, str, str, str], bool] = {}
    # iterate over donors from earliest submission to latest, later overriding earlier
    for data_node_id, submitter_id, index_local_case_id, donor in prior_donors:
        pseudonym = index_local_case_id if donor.relation == Relation.index_ else donor.pseudonym
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
    quarter_start_date, quarter_end_date = quarter_date_bounds(year=year, quarter=quarter)

    node_submitter_id_combos: set[tuple[str | None, str | None]] = set()
    with database.transaction() as session:
        # number_of_end-to-end_tests
        stmt_number_of_end_to_end_tests = (
            select(Submission.data_node_id, Submission.submitter_id, sqlfn.count(1))
            .where(Submission.submission_uploaded_date.between(quarter_start_date, quarter_end_date))  # type: ignore[union-attr]
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
            .where(Submission.submission_uploaded_date.between(quarter_start_date, quarter_end_date))  # type: ignore[union-attr]
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
            .where(Submission.submission_uploaded_date.between(quarter_start_date, quarter_end_date))  # type: ignore[union-attr]
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
            .where(Submission.submission_uploaded_date.between(quarter_start_date, quarter_end_date))  # type: ignore[union-attr]
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
            .where(Submission.submission_uploaded_date.between(quarter_start_date, quarter_end_date))  # type: ignore[union-attr]
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
        # skip over combos with no recorded submitter or data node ID
        for data_node_id, submitter_id in sorted(
            filter(lambda p: all(map(lambda i: i is not None, p)), node_submitter_id_combos)
        ):
            if (data_node_id is None) or (submitter_id is None):
                # needed for type checking
                raise ValueError("Encountered null data node ID or submitter ID after they should have been filtered")

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


def _detailed_qc_passed_state(submission: Submission) -> DetailedQCPassedReportState:
    """Report the stored detailed QC verdict, which is unset until detailed QC has run."""
    if submission.detailed_qc_passed is None:
        return DetailedQCPassedReportState.NOT_PERFORMED
    return DetailedQCPassedReportState.YES if submission.detailed_qc_passed else DetailedQCPassedReportState.NO


def _dump_dataset_report(
    output_path: Path, database: SubmissionDb, year: int, quarter: int, with_submission_ids: bool = False
) -> None:
    quarter_start_date, quarter_end_date = quarter_date_bounds(year=year, quarter=quarter)

    with database.transaction() as session:
        query_quarter_submissions = (
            select(Submission).where(Submission.submission_uploaded_date.between(quarter_start_date, quarter_end_date))  # type: ignore[union-attr]
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
            + (["submission_id"] if with_submission_ids else [])
        )
        for submission in submissions:
            detailed_qc_passed = _detailed_qc_passed_state(submission)
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
                + ([submission.id] if with_submission_ids else [])
            )


def _dump_qc_report(
    output_path: Path, database: SubmissionDb, year: int, quarter: int, with_submission_ids: bool = False
) -> None:
    quarter_start_date, quarter_end_date = quarter_date_bounds(year=year, quarter=quarter)

    with database.transaction() as session:
        # The quarterly QC report lists submissions that failed the detailed QC (a computed
        # metric below the required threshold) as well as submissions that passed but carry a
        # deviation of more than PCT_DEV_CUTOFF percent between the computed and provided values
        # on any metric. BfArM requires both to be reported; the latter does not count as a
        # failure. The detailed_qc_passed and deviation_exceeds_tolerance columns tell the two
        # reasons apart per row.
        submission_ids_with_deviation = select(DetailedQCResult.submission_id).where(
            sa.or_(
                sa.not_(DetailedQCResult.percent_bases_above_quality_threshold_passed_qc),  # type: ignore[call-overload]
                sa.not_(DetailedQCResult.mean_depth_of_coverage_passed_qc),  # type: ignore[call-overload]
                sa.not_(DetailedQCResult.targeted_regions_above_min_coverage_passed_qc),  # type: ignore[call-overload]
            )
        )
        query_submissions_to_report = (
            select(Submission)
            .where(Submission.submission_uploaded_date.between(quarter_start_date, quarter_end_date))  # type: ignore[union-attr]
            .filter(
                sa.or_(
                    sa.not_(Submission.detailed_qc_passed),  # type: ignore[call-overload]
                    Submission.id.in_(submission_ids_with_deviation),  # type: ignore[attr-defined]
                )
            )
        )
        submissions_to_report = session.exec(query_submissions_to_report).all()

        subquery_submissions_to_report = query_submissions_to_report.subquery()
        query_reports_to_include = select(DetailedQCResult).join(
            subquery_submissions_to_report,
            subquery_submissions_to_report.c.id == DetailedQCResult.submission_id,
        )
        reports_to_include = session.exec(query_reports_to_include).all()

        subquery_reports_to_include = query_reports_to_include.subquery()
        query_relevant_donors = select(Donor).join(
            subquery_reports_to_include,
            (subquery_reports_to_include.c.submission_id == Donor.submission_id)  # type: ignore[arg-type]
            & (subquery_reports_to_include.c.pseudonym == Donor.pseudonym),
        )
        relevant_donors = session.exec(query_relevant_donors).all()

    # Note: for future, this join only works because *both* DetailedQCResult
    # and Donor store the index donor's pseudonym as "index"; the identifier the
    # submitter uses for the index donor is submissions.local_case_id
    idpseudo2relation = {(donor.submission_id, donor.pseudonym): donor.relation for donor in relevant_donors}
    id2submission = {submission.id: submission for submission in submissions_to_report}
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
                "detailed_qc_passed",
                "deviation_exceeds_tolerance",
            ]
            + (["submission_id"] if with_submission_ids else [])
        )

        for report in reports_to_include:
            submission = id2submission[report.submission_id]
            relation = Relation(idpseudo2relation[(report.submission_id, report.pseudonym)])
            deviation_exceeds_tolerance = any(
                abs(deviation) > PCT_DEV_CUTOFF
                for deviation in (
                    report.percent_bases_above_quality_threshold_percent_deviation,
                    report.mean_depth_of_coverage_percent_deviation,
                    report.targeted_regions_above_min_coverage_percent_deviation,
                )
            )
            writer.writerow(
                [
                    submission.data_node_id,
                    quarter,
                    year,
                    submission.submitter_id,
                    submission.submission_uploaded_date,
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
                    _detailed_qc_passed_state(submission),
                    "yes" if deviation_exceeds_tolerance else "no",
                ]
                + ([report.submission_id] if with_submission_ids else [])
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
    help="Directory to output TSV files. Defaults to current directory.",
)
@click.option(
    "--with-submission-ids",
    "with_submission_ids",
    is_flag=True,
    default=False,
    help="Append a 'submission_id' column to per-submission tables (dataset and detailed QC).",
)
@click.pass_context
def quarterly(
    ctx: click.Context, year: int | None, quarter: int | None, output_directory: Path | None, with_submission_ids: bool
):
    """
    Generate the tables for the quarterly report.
    """
    db = ctx.obj["db_url"]
    grz_id = ctx.obj["grz_id"]
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

    overview_output_path = output_directory / f"1-Gesamtübersicht_{grz_id}_{quarter}_{year}.tsv"
    _dump_overview_report(overview_output_path, submission_db, year, quarter)

    dataset_output_path = output_directory / f"2-Infos_zu_Datensätzen_{grz_id}_{quarter}_{year}.tsv"
    _dump_dataset_report(dataset_output_path, submission_db, year, quarter, with_submission_ids=with_submission_ids)

    qc_output_path = output_directory / f"3-Detailprüfung_{grz_id}_{quarter}_{year}.tsv"
    _dump_qc_report(qc_output_path, submission_db, year, quarter, with_submission_ids=with_submission_ids)


def _utc(timestamp: datetime.datetime) -> datetime.datetime:
    """SQLite returns naive timestamps; they are stored in UTC."""
    return timestamp if timestamp.tzinfo else timestamp.replace(tzinfo=datetime.UTC)


def _monitoring_payload(
    grz_id: str, month: str, le_ids: typing.Iterable[str], submissions: typing.Iterable[Submission]
) -> dict[str, typing.Any]:
    """Cumulative per-LE numbers for the GRZ monitoring sheet, as of the end of *month* (``YYYY-MM``).

    - ``cases``: cases whose initial submission passed basic QC and had its Prüfbericht submitted (state
      Reported) by then. A case has at most one such submission, so this counts cases, not submissions:
      follow-up, addition and correction submissions belong to a case already counted, and a rejected
      initial submission is no case. The QC flag matters because a Prüfbericht can be reported for a
      rejected submission too: older databases hold negative Prüfberichte with a Reported state.
    - ``avg_gb``: mean size of those initial submissions, in GB (10^9 bytes).
    - ``active``: 1 if the LE uploaded a non-test submission during *month*.
    - ``tb_consented`` / ``tb_unconsented``: archived, not deleted submissions by research consent, in TB (10^12
      bytes). A missing consent counts as unconsented, as in ``grzctl archive``.

    Every LE in *le_ids* gets a row, also with zero cases, so the sheet sees the month as complete.
    """
    first = datetime.date.fromisoformat(f"{month}-01")
    after = datetime.date(first.year + first.month // 12, first.month % 12 + 1, 1)
    cutoff = datetime.datetime.combine(after, datetime.time(), datetime.UTC)

    def reached(log: typing.Iterable[tuple[str, datetime.datetime]], kind: str) -> bool:
        return any(entry == kind and _utc(timestamp) < cutoff for entry, timestamp in log)

    cases: defaultdict[str, int] = defaultdict(int)
    sizes: defaultdict[str, list[int]] = defaultdict(list)
    active: set[str] = set()
    archived_bytes = {True: 0, False: 0}
    for submission in submissions:
        le_id = submission.submitter_id
        if le_id is None or submission.submission_type == SubmissionType.test:
            continue
        states = [(state.state, state.timestamp) for state in submission.states]
        if (
            submission.submission_type == SubmissionType.initial
            and submission.basic_qc_passed
            and reached(states, SubmissionStateEnum.REPORTED)
        ):
            cases[le_id] += 1
            if submission.submission_size is not None:
                sizes[le_id].append(submission.submission_size)
        if submission.submission_uploaded_date and first <= submission.submission_uploaded_date < after:
            active.add(le_id)
        changes = [(change.change, change.timestamp) for change in submission.changes]
        deleted = reached(changes, ChangeRequestEnum.DELETE)
        if reached(states, SubmissionStateEnum.ARCHIVED) and not deleted:
            archived_bytes[bool(submission.consented)] += submission.submission_size or 0

    rows = []
    for le_id in sorted(set(le_ids) | set(cases) | active):
        row: dict[str, typing.Any] = {"le_id": le_id, "cases": cases[le_id], "active": int(le_id in active)}
        if sizes[le_id]:
            row["avg_gb"] = sum(sizes[le_id]) / len(sizes[le_id]) / 1e9
        rows.append(row)
    return {
        "grz": grz_id,
        "month": month,
        "reported_by": "grzctl report monitoring",
        "tb_consented": archived_bytes[True] / 1e12,
        "tb_unconsented": archived_bytes[False] / 1e12,
        "rows": rows,
    }


@report.command()
@click.option(
    "--month",
    default=lambda: datetime.date.today().strftime("%Y-%m"),
    show_default="current month",
    help="Report the state at the end of this month (YYYY-MM).",
)
@click.option(
    "--url",
    envvar="GRZ_MONITORING_URL",
    help="Intake URL of the GRZ monitoring sheet (env GRZ_MONITORING_URL). Without it, the report is only printed.",
)
@click.option(
    "--token",
    envvar="GRZ_MONITORING_TOKEN",
    help="This GRZ's intake token. Prefer the environment variable GRZ_MONITORING_TOKEN.",
)
@click.option("--dry-run", is_flag=True, help="Let the intake check the report without storing it.")
@click.pass_context
def monitoring(ctx: click.Context, month: str, url: str | None, token: str | None, dry_run: bool):
    """
    Send cumulative per-LE numbers to the GRZ monitoring sheet.

    Meant for a daily or weekly cron job. Sending a month again replaces its numbers in the sheet.
    """
    if not re.fullmatch(r"\d{4}-(0[1-9]|1[0-2])", month):
        raise click.BadParameter("must be YYYY-MM", param_hint="--month")

    submission_db = get_submission_db_instance(ctx.obj["db_url"])
    with submission_db.transaction() as session:
        submissions = session.exec(
            select(Submission).options(selectinload(Submission.states), selectinload(Submission.changes))  # type: ignore[arg-type]
        ).all()
        payload = _monitoring_payload(ctx.obj["grz_id"], month, ctx.obj["le_ids"], submissions)

    if url is None:
        click.echo(json.dumps(payload, indent=2))
        return
    if not token:
        raise click.UsageError("--url needs --token or GRZ_MONITORING_TOKEN.")

    # Apps Script answers with a redirect; requests follows it as GET, which is what the intake expects.
    response = requests.post(url, json={**payload, "token": token, "dry_run": dry_run}, timeout=60)
    response.raise_for_status()
    result = response.json()
    click.echo(json.dumps(result))
    if result.get("status") != 200:
        raise click.ClickException(f"intake rejected the report: {result.get('error', result)}")
