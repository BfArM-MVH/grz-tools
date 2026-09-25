"""Inbox resolution for grzctl commands.

A submission is uploaded to the S3 inbox(es) of one submitter (LE).
Commands that read or clean a single submission resolve the inbox from, in order:
an explicit ``--inbox``, the inbox recorded in the database for that submission, the
submitter's only inbox, and, when asked, by scanning the submitter's inboxes for the submission.
"""

import click
from grz_common.workers.download import query_submissions
from grz_db.models.submission import SubmissionDb

from ..models.config import GrzctlConfig


def db_inbox(db_service: SubmissionDb | None, submission_id: str) -> str | None:
    """Return the inbox the database recorded for ``submission_id``, if any.

    A database that is not available, or does not answer, yields ``None`` rather than an error.
    The inbox simply is not known, and the remaining sources still resolve.
    """
    if db_service is None:
        return None
    try:
        submission = db_service.get_submission(submission_id)
    except Exception:
        return None
    return submission.inbox if submission is not None and submission.inbox else None


def _only_inbox(configuration: GrzctlConfig, submitter_id: str) -> str | None:
    """Return the submitter's only inbox, or ``None`` when it has none or several."""
    entry = configuration.leistungserbringer.get(submitter_id)
    if entry is None or len(entry.inbox_buckets) != 1:
        return None
    return next(iter(entry.inbox_buckets))


def _scan_inbox(configuration: GrzctlConfig, submitter_id: str, submission_id: str) -> str | None:
    """Return the sole inbox of ``submitter_id`` that still holds ``submission_id``.

    Cleaning keeps a (redacted) ``metadata.json`` marker in the inbox.
    A cleaned submission is therefore still found here.
    Found in several inboxes, or in none, yields ``None``.
    Both cases are ambiguous and need an explicit ``--inbox``.
    """
    entry = configuration.leistungserbringer.get(submitter_id)
    if entry is None:
        return None
    matches = [
        inbox_name
        for inbox_name in entry.inbox_buckets
        if submission_id in _inbox_listing(configuration, submitter_id, inbox_name)
    ]
    return matches[0] if len(matches) == 1 else None


def resolve_inbox(  # noqa: PLR0913
    configuration: GrzctlConfig,
    *,
    submitter_id: str,
    submission_id: str | None = None,
    inbox_name: str | None = None,
    db_service: SubmissionDb | None = None,
    scan: bool = False,
) -> str | None:
    """Resolve the inbox a submission-scoped command should operate on.

    Resolution order:

    1. ``inbox_name``, when given.
    2. The inbox the database recorded for ``submission_id``, when ``db_service`` is set.
    3. The submitter's only inbox.
    4. With ``scan``, the sole inbox still holding ``submission_id`` in S3.

    The bucket always follows from the inbox via :meth:`GrzctlConfig.inbox_target`.

    :param configuration: The grzctl configuration, which names the submitter's inboxes.
    :param submitter_id: The submitter (LE) the submission belongs to.
    :param submission_id: The submission, for the database and scan lookups.
    :param inbox_name: An explicitly given inbox; it wins over every other source.
    :param db_service: Submission database to look the recorded inbox up in.
    :param scan: Whether to scan the submitter's inboxes for ``submission_id`` when nothing
        else resolved.
    :returns: The inbox name, or ``None`` when nothing resolves.
    """
    if inbox_name is not None:
        return inbox_name
    if db_service is not None and submission_id is not None:
        recorded = db_inbox(db_service, submission_id)
        if recorded is not None:
            return recorded
    only = _only_inbox(configuration, submitter_id)
    if only is not None:
        return only
    if scan and submission_id is not None:
        return _scan_inbox(configuration, submitter_id, submission_id)
    return None


def require_inbox(  # noqa: PLR0913
    configuration: GrzctlConfig,
    *,
    submitter_id: str,
    submission_id: str | None = None,
    inbox_name: str | None = None,
    db_service: SubmissionDb | None = None,
    scan: bool = False,
    hint: str = "Pass --inbox.",
    exc_type: type[Exception] = click.ClickException,
) -> str:
    """Resolve the inbox a command should operate on, or abort.

    Resolution order is that of :func:`resolve_inbox`.
    When nothing resolves, it raises an error that names what is known about the
    submitter's inboxes and repeats *hint*, which names the alternatives that resolve.
    Commands scoped to a single submission pass its ``submission_id``.

    :param configuration: The grzctl configuration, which names the submitter's inboxes.
    :param submitter_id: The submitter (LE) the submission belongs to.
    :param submission_id: The submission, for the database and scan lookups.
    :param inbox_name: An explicitly given inbox; it wins over every other source.
    :param db_service: Submission database to look the recorded inbox up in.
    :param scan: Whether to scan the submitter's inboxes for ``submission_id`` when nothing else resolved.
    :param hint: Advice appended to the error, naming options that resolve.
    :param exc_type: The error to raise; ``click.UsageError`` for submitter-scoped commands.
        ``decrypt`` passes :class:`~grz_common.exceptions.ConfigurationError`, which its ``DbContext`` records.
    :returns: The resolved inbox name.
    :raises click.ClickException: if no inbox resolves, or *exc_type* if given.
    """
    resolved = resolve_inbox(
        configuration,
        submitter_id=submitter_id,
        submission_id=submission_id,
        inbox_name=inbox_name,
        db_service=db_service,
        scan=scan,
    )
    if resolved is not None:
        return resolved
    entry = configuration.leistungserbringer.get(submitter_id)
    if submission_id is not None:
        where = f"for submission {submission_id} of submitter {submitter_id}"
    else:
        where = f"for submitter {submitter_id}"
    if entry is None:
        detail = "the submitter has no inbox in the configuration"
    elif len(entry.inbox_buckets) == 1:
        detail = "the submission is not present in the submitter's only inbox"
    else:
        detail = (
            f"the submitter has several inboxes ({', '.join(sorted(entry.inbox_buckets))}), and none is unambiguous"
        )
    raise exc_type(f"Cannot resolve an inbox {where}: {detail}. {hint}")


def scan_inbox(configuration: GrzctlConfig, submitter_id: str, submission_id: str) -> str | None:
    """Return the sole inbox of ``submitter_id`` that still holds ``submission_id``, if unambiguous.

    Unlike :func:`resolve_inbox` this only consults S3.
    It is for filling a recorded inbox where there is none yet, e.g. during backfill.
    Listings are cached per inbox, so backfilling many submissions lists each inbox once.
    """
    return _scan_inbox(configuration, submitter_id, submission_id)


def _inbox_listing(configuration: GrzctlConfig, submitter_id: str, inbox_name: str) -> frozenset[str]:
    """Submission IDs found in one inbox, listing each inbox at most once per process."""
    key = (submitter_id, inbox_name)
    if key not in _inbox_listing_cache:
        s3_options = configuration.inbox_target(submitter_id=submitter_id, inbox_name=inbox_name).s3
        _inbox_listing_cache[key] = frozenset(
            summary.submission_id for summary in query_submissions(s3_options, show_cleaned=True)
        )
    return _inbox_listing_cache[key]


_inbox_listing_cache: dict[tuple[str, str], frozenset[str]] = {}
