# Case tracking

A _case_ groups the submissions belonging to one patient. Every non-`test` submission is
linked to a case during `populate` or `backfill`, and the link is stored on
`submissions.case_id`. This page covers what a case is, what the upgrade does to existing
data, and how to repair a link by hand.

> ⚠️ **The upgrade reports unlinkable data once, to its log.** `grzctl db upgrade` names the
> submitter-local case keys it refused to group, and nothing else records them afterwards.
> Capture that output. [Upgrading an existing database](#upgrading-an-existing-database)
> explains what to do with it, and how to re-derive the list in SQL if it was lost.

- [Concepts](#concepts)
- [One initial submission per case](#one-initial-submission-per-case)
- [Upgrading an existing database](#upgrading-an-existing-database)
- [Repairing links by hand](#repairing-links-by-hand)
- [Populate and backfill](#populate-and-backfill)
- [Command reference](#command-reference)

## Concepts

A case has three identifying columns, and they are not equally authoritative:

| Column          | Description                                                          |
| --------------- | -------------------------------------------------------------------- |
| `psn`           | The RKI pseudonym. The authoritative identity, unique once assigned. |
| `submitter_id`  | The Leistungserbringer that sent the submissions.                    |
| `local_case_id` | The case identifier that the submitter uses for the patient.         |

`psn` is the identity. `submitter_id` and `local_case_id` are a _resolution key_, used to
locate a case before a `psn` exists. A partial unique index keeps the pair unique wherever
both halves are present. Neither half is required, so a future flow can resolve a case by
`psn` alone.

Two naming traps are worth knowing before you read any query:

- **`submissions.pseudonym` is not the `psn`.** That column holds the submitter's
  `localCaseId`. The quarterly report also uses it to stand in for the index donor, whose
  own `donorPseudonym` is the literal `"index"`.
- **A redaction placeholder is not a case key.** Archived documents carry either an empty
  `localCaseId` or `REDACTED_LOCAL_CASE_ID`. Neither identifies a patient, so neither opens
  or resolves a case. Restore the submitter's value first when you work from an archived
  copy.

`test` submissions are never case-tracked. Their `case_id` stays NULL permanently, and the
`db case` commands refuse to link them.

## One initial submission per case

A case may have at most one `initial` submission that passed basic QC. A partial unique
index enforces it, so the rule holds no matter which door a write arrives at.

Linking itself is permissive. Several `initial` submissions may share a case while their
basic QC is still pending, because the data alone cannot tell a re-upload from a genuine
duplicate. The rejection lands when a second one tries to _pass_.

`grzctl validate` asks the question before doing the work. When the case already has a
QC-passed initial, validation stops, the submission's `basic_qc_passed` is set to false, and
the state log records an `ERROR` with failure reason `duplicate_initial`. The same handling
applies to a rival that passes while this submission is being validated.

With `--no-update-db`, nothing is written, so the check only logs a warning.

## Upgrading an existing database

The cases migration adds the `cases` table and `submissions.case_id`, groups existing
submissions by `(submitter_id, pseudonym)`, and creates one case per group.

**Run the upgrade with its output captured:**

```bash
grzctl --config "$CONFIG" db upgrade 2>&1 | tee db-upgrade-$(date +%F).log
```

`--config` is an option of `grzctl` itself, so it comes before `db`. Every example below
omits it for brevity.

Some rows are deliberately left unlinked: rows missing either half of the key, `test`
submissions, and rows whose key is shared by more than one QC-passed `initial` submission.
That last group is the one that needs a human. A key naming two QC-passed initials names two
patients, so no case is created for it and every submission carrying it keeps `case_id` NULL:

```
WARNING  [alembic.runtime.migration] cases backfill: 2 submitter-local case key(s) are
shared by more than one QC-passed 'initial' submission, so they identify no single patient.
No case is created for them and every submission carrying one keeps case_id NULL:
  (260914050, patient-0042): 2 QC-passed initial submissions
  (260914050, patient-0117): 3 QC-passed initial submissions
```

NULL is the "not yet resolved" state here, not a failure. Nothing is lost, and a later flow
can still resolve these submissions by `psn`. Merging two patients into one case would have
been the destructive choice, because unpicking it means one `db case unlink` per submission
and the key that caused the merge tells you nothing about which submission belongs where.

If the output was not captured, the same list can be re-derived:

```sql
SELECT submitter_id, pseudonym, count(*)
FROM submissions
WHERE submitter_id IS NOT NULL
  AND pseudonym IS NOT NULL
  AND pseudonym NOT IN ('', 'REDACTED_LOCAL_CASE_ID')
  AND submission_type = 'initial'
  AND basic_qc_passed IS TRUE
GROUP BY submitter_id, pseudonym
HAVING count(*) > 1;
```

Resolving one of these means deciding which submissions belong to which patient, then
creating a case per patient and linking its submissions:

```bash
grzctl --config "$CONFIG" db case create 260914050 patient-0042 --psn RKI-000123
grzctl --config "$CONFIG" db case relink 260914050_2025-03-11_a1b2c3d4 7
```

## Repairing links by hand

Cases are normally created by `populate`. The commands below exist for repair.

```bash
db case list                  # every case with its linked-submission count
db case show 7                # one case and its submissions, with basic QC state
db case show 7 --json         # same, machine-readable
```

`db case show` lists each submission's basic QC state, which is how you see which `initial`
holds the case's QC-passed slot.

```bash
db case create SUBMITTER_ID LOCAL_CASE_ID [--psn PSN]
db case modify 7 psn RKI-000123
db case relink SUBMISSION_ID CASE_ID
db case unlink SUBMISSION_ID
db case delete 7
```

Four behaviours worth knowing before you use them:

- `relink` moves a submission between cases. It cannot leave it linked to none, and it still
  refuses a move that would give the target case a second QC-passed initial.
- `unlink` is the only route back to NULL.
- Either one leaves the vacated case as it is, so a case may end up empty.
- `delete` refuses a case that still has submissions linked. Unlink them first.

## Populate and backfill

`populate` and `backfill` resolve the case from the metadata and write the link like any
other missing value. Two rules differ from ordinary columns:

- `--allow-overwrite` cannot name the case link. Replacing an existing link would undo a
  deliberate `db case relink`, so it holds the whole submission back until `--force` is
  passed.
- `--ignore-field case_id` skips case resolution altogether.

A case that cannot be resolved does not cost the submission the rest of its metadata.
Everything else is written, and `backfill` reports the submission under its own counter:

```
Done. Updated: 1834
  Up to date: 12
  Not in bucket (split consent): 3
  Would overwrite (needs --force): 0
  Case link unresolved: 2 (also counted above)
  Errors: 0
```

Re-running links those submissions once the cause is cleared.

Note that `download --populate` does not surface an unresolvable case. The download
succeeds, the submission is left unlinked, and only a later `backfill` or `populate` reports
it.

## Command reference

| Command                                                 | Purpose                                          |
| ------------------------------------------------------- | ------------------------------------------------ |
| `db case list [--json]`                                 | All cases with their linked-submission counts.   |
| `db case show CASE_ID [--json]`                         | One case, its columns and its submissions.       |
| `db case create SUBMITTER_ID LOCAL_CASE_ID [--psn PSN]` | Create a case by hand.                           |
| `db case modify CASE_ID KEY VALUE`                      | Change `local_case_id`, `psn` or `submitter_id`. |
| `db case relink SUBMISSION_ID CASE_ID`                  | Move a submission to a different case.           |
| `db case unlink SUBMISSION_ID`                          | Detach a submission from its case.               |
| `db case delete CASE_ID`                                | Delete a case that has no submissions linked.    |
