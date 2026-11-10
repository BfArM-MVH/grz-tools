"""Command for generating various reports related to GRZ activities."""

import csv
import datetime
import itertools
import json
import logging
import statistics
import typing
from collections import defaultdict
from enum import StrEnum
from operator import attrgetter
from pathlib import Path

import click
import requests
import sqlalchemy as sa
from grz_common.transfer import UploadTiming
from grz_db.models.submission import (
    ChangeRequestEnum,
    ChangeRequestLog,
    DetailedQCResult,
    Donor,
    Submission,
    SubmissionDb,
    SubmissionStateEnum,
    SubmissionStateLog,
)
from grz_pydantic_models.dates import date_to_quarter_year, month_date_bounds, quarter_date_bounds
from grz_pydantic_models.submission.metadata import GenomicStudyType, Relation, SubmissionType
from grz_pydantic_models.submission.thresholds import PCT_DEV_CUTOFF
from sqlalchemy import func as sqlfn
from sqlmodel import Session, select

from .. import __version__
from ..commands import grzctl_configuration
from ..models.config import GrzctlConfig
from .db.cli import get_submission_db_instance

log = logging.getLogger(__name__)


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
        "submitter_ids": list(configuration.leistungserbringer),
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


def _monitoring_payload(
    session: Session, grz_id: str, month: datetime.date, submitter_ids: typing.Iterable[str]
) -> dict[str, typing.Any]:
    """Per-LE numbers for the GRZ monitoring sheet, as of the end of the month of *month*.

    Every LE in *submitter_ids* gets a row, also without submissions, so the sheet sees the month as
    complete, and every row has the same keys: a count that has nothing to count is ``0``, a statistic
    without values ``None``. Sizes are in bytes; the sheet scales them. Months end at midnight UTC, as
    the database timestamps are in UTC. Only the state and change-request logs are read as of the
    month's end; QC result, size and consent are the submissions' current values.

    Cumulative, so the sheet derives monthly increments from them, and something recorded late shows up
    in a later increment instead of getting lost:

    - ``submissions``: non-test submissions uploaded by then. An LE was active in a month if this grew.
    - ``cases``: initial submissions that passed basic QC and had their Prüfbericht submitted (state
      Reported) by then. An unset QC flag counts if the submission was archived, since ``grzctl validate``
      has only set it from 2026-07 on and a failed submission is not archived. Follow-up, addition and
      correction submissions belong to a case already counted, and a rejected initial submission is no
      case. The database allows a case only one initial submission with the flag set, so this counts
      cases, not submissions; for an unset flag, which that index does not cover, it assumes that a case
      has no second archived initial submission.
    - ``cases_with_size`` and ``case_bytes``: how many of those initial submissions have a known size,
      and their total size. The sheet divides the one by the other.
    - ``archived`` and ``archived_with_size``: the LE's archived submissions (state Archived, or
      Processed from ``grzctl process``) that were neither deleted nor transferred to another GRZ, and
      how many of them have a known size. A missing size is no size of 0, so the sheet can tell the
      volume below is incomplete.
    - ``archived_bytes_consented``, ``archived_bytes_non_consented`` and
      ``archived_bytes_consent_unknown``: their total known size by research consent. Metadata without
      a research consent is stored as ``False``, so it is non-consented, as in ``grzctl archive``.
      Unknown means the database holds no value.
    - ``timed_uploads``, ``timed_upload_bytes`` and ``timed_upload_seconds``: number, total size and
      total duration of the LE's non-test uploads with a known duration that finished by then, as
      ``grzctl download`` or ``grzctl process`` recorded them (:attr:`Submission.upload_started_at`).
      A resumed upload and an upload without the uploader's ``version`` marker have none: the first
      took longer than its last run, and for the second a resumption cannot be detected
      (``upload_resumed`` is NULL). The sheet divides the increments of the bytes by those of the
      seconds to get the mean upload speed.

    Of the month alone, because they cannot be accumulated:

    - ``upload_speed_{median,min,max,sd}_mb_per_s``: median, minimum, maximum and population standard
      deviation of the upload speed of those timed uploads that finished during the month.
    """
    first, last = month_date_bounds(month.year, month.month)
    start = datetime.datetime.combine(first, datetime.time(), datetime.UTC)
    cutoff = datetime.datetime.combine(last + datetime.timedelta(days=1), datetime.time(), datetime.UTC)
    real = sa.and_(
        Submission.submitter_id.is_not(None),  # type: ignore[union-attr]
        Submission.submission_type != SubmissionType.test,  # type: ignore[arg-type]
    )

    def logged(log: typing.Any, column: typing.Any, kind: str) -> sa.ColumnElement[bool]:
        """A log entry of *kind* for the submission before the cutoff."""
        return sa.exists().where(log.submission_id == Submission.id, column == kind, log.timestamp < cutoff)

    reported = logged(SubmissionStateLog, SubmissionStateLog.state, SubmissionStateEnum.REPORTED)
    # ``grzctl process`` records Processed once it has archived a submission, instead of Archived
    archived = sa.or_(
        logged(SubmissionStateLog, SubmissionStateLog.state, SubmissionStateEnum.ARCHIVED),
        logged(SubmissionStateLog, SubmissionStateLog.state, SubmissionStateEnum.PROCESSED),
    )
    # a deleted submission is gone, a transferred one is held by another GRZ: neither is archived here anymore
    gone = sa.or_(
        logged(ChangeRequestLog, ChangeRequestLog.change, ChangeRequestEnum.DELETE),
        logged(ChangeRequestLog, ChangeRequestLog.change, ChangeRequestEnum.TRANSFER),
    )
    with_size = sqlfn.count(Submission.submission_size)  # type: ignore[arg-type]
    total_size = sqlfn.coalesce(sqlfn.sum(Submission.submission_size), 0)

    submissions: dict[str, int] = dict(
        session.execute(
            select(Submission.submitter_id, sqlfn.count(1))
            .where(real, Submission.submission_uploaded_date <= last)  # type: ignore[operator]
            .group_by(Submission.submitter_id)
        )
        .tuples()
        .all()
    )
    cases = {
        submitter_id: (count, sized, int(total))
        for submitter_id, count, sized, total in session.execute(
            select(Submission.submitter_id, sqlfn.count(1), with_size, total_size)
            .where(
                real,
                Submission.submission_type == SubmissionType.initial,  # type: ignore[arg-type]
                # an unset flag (set by ``grzctl validate`` only from 2026-07 on) counts once archived:
                # a submission that failed QC is not archived
                sa.or_(Submission.basic_qc_passed.is_(True), sa.and_(Submission.basic_qc_passed.is_(None), archived)),  # type: ignore[union-attr]
                reported,
            )
            .group_by(Submission.submitter_id)
        ).all()
    }
    archived_count: defaultdict[str, int] = defaultdict(int)
    archived_sized: defaultdict[str, int] = defaultdict(int)
    archived_bytes: defaultdict[tuple[str, bool | None], int] = defaultdict(int)
    for submitter_id, consented, count, sized, total in session.execute(
        select(Submission.submitter_id, Submission.consented, sqlfn.count(1), with_size, total_size)  # type: ignore[call-overload]
        .where(real, archived, sa.not_(gone))
        .group_by(Submission.submitter_id, Submission.consented)  # type: ignore[arg-type]
    ).all():
        archived_count[submitter_id] += count
        archived_sized[submitter_id] += sized
        archived_bytes[submitter_id, consented] = int(total)

    timed_uploads: defaultdict[str, int] = defaultdict(int)
    timed_bytes: defaultdict[str, int] = defaultdict(int)
    timed_seconds: defaultdict[str, float] = defaultdict(float)
    speeds: defaultdict[str, list[float]] = defaultdict(list)
    for submitter_id, size, started, finished, in_month in session.execute(
        select(  # type: ignore[call-overload]
            Submission.submitter_id,
            Submission.submission_size,
            Submission.upload_started_at,
            Submission.upload_finished_at,
            Submission.upload_finished_at >= start,  # type: ignore[operator]
        ).where(
            real,
            Submission.submission_size > 0,  # type: ignore[operator]
            Submission.upload_resumed.is_(False),  # type: ignore[union-attr]
            Submission.upload_started_at.is_not(None),  # type: ignore[union-attr]
            Submission.upload_finished_at < cutoff,  # type: ignore[operator]
        )
    ).all():
        seconds = UploadTiming(started, finished, upload_resumed=False).duration.total_seconds()
        if seconds > 0:
            timed_uploads[submitter_id] += 1
            timed_bytes[submitter_id] += size
            timed_seconds[submitter_id] += seconds
            if in_month:
                speeds[submitter_id].append(size / 1e6 / seconds)

    statistic: dict[str, typing.Callable[[list[float]], float]] = {
        "median": statistics.median,
        "min": min,
        "max": max,
        "sd": statistics.pstdev,
    }
    rows = []
    for submitter_id in sorted(
        set(submitter_ids) | set(submissions) | set(cases) | set(archived_count) | set(timed_uploads)
    ):
        count, sized, total = cases.get(submitter_id, (0, 0, 0))
        mb_per_s = speeds.get(submitter_id, [])
        rows.append(
            {
                "submitter_id": submitter_id,
                "submissions": submissions.get(submitter_id, 0),
                "cases": count,
                "cases_with_size": sized,
                "case_bytes": total,
                "archived": archived_count[submitter_id],
                "archived_with_size": archived_sized[submitter_id],
                "archived_bytes_consented": archived_bytes[submitter_id, True],
                "archived_bytes_non_consented": archived_bytes[submitter_id, False],
                "archived_bytes_consent_unknown": archived_bytes[submitter_id, None],
                "timed_uploads": timed_uploads[submitter_id],
                "timed_upload_bytes": timed_bytes[submitter_id],
                "timed_upload_seconds": round(timed_seconds[submitter_id]),
                **{
                    f"upload_speed_{name}_mb_per_s": round(f(mb_per_s), 2) if mb_per_s else None
                    for name, f in statistic.items()
                },
            }
        )
    return {
        "grz": grz_id,
        "month": month.strftime("%Y-%m"),
        "reported_by": f"grzctl {__version__} report monitoring",
        "rows": rows,
    }


@report.command()
@click.option(
    "--month",
    type=click.DateTime(formats=["%Y-%m"]),
    show_default="the running month and the one before",
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
def monitoring(ctx: click.Context, month: datetime.datetime | None, url: str | None, token: str | None, dry_run: bool):
    """
    Send per-LE numbers to the GRZ monitoring sheet, one report per month.

    Meant for a daily cron job. Sending a month again replaces its numbers in the sheet. Without
    --month, each run sends the running month and the one before: the month before gets its final
    figures, and its upload speed statistics still collect uploads downloaded late.
    """
    today = datetime.datetime.now(datetime.UTC).date()  # the month ends of the report are in UTC
    if month and month.date() > today:
        raise click.BadParameter(f"{month:%Y-%m} has not begun yet.", param_hint="--month")
    months = [month.date()] if month else [today.replace(day=1) - datetime.timedelta(days=1), today]
    submission_db = get_submission_db_instance(ctx.obj["db_url"])
    with submission_db.transaction() as session:
        payloads = [_monitoring_payload(session, ctx.obj["grz_id"], m, ctx.obj["submitter_ids"]) for m in months]

    if url is None:
        click.echo(json.dumps(payloads, indent=2))
        return
    if not token:
        raise click.UsageError("--url needs --token or GRZ_MONITORING_TOKEN.")

    for payload in payloads:
        # Apps Script answers with a redirect; requests follows it as GET, which is what the intake expects.
        try:
            response = requests.post(url, json={**payload, "token": token, "dry_run": dry_run}, timeout=60)
            response.raise_for_status()
            result = response.json()
        except requests.RequestException as e:
            raise click.ClickException(f"intake not reachable: {e}") from e
        click.echo(json.dumps(result))
        if not isinstance(result, dict):
            raise click.ClickException(f"intake answered the report of {payload['month']} with no JSON object")
        if result.get("status") != 200:
            raise click.ClickException(
                f"intake rejected the report of {payload['month']}: {result.get('error', result)}"
            )
