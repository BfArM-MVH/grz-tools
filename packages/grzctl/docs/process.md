# Processing a submission

`grzctl process` runs a submission through the whole pipeline in a single streaming
pass: download the metadata, validate and re-encrypt each file, stage the result in
the _interrogation bucket_, decide whether the submission goes through detailed QC,
then commit everything to the target archive bucket and clean up the inbox. It also
drives the submission's database state and its QC selection.

The one thing to know as an operator: the pipeline guesses the detailed-QC decision
early, before the real decision exists, and prefetches a decrypted copy of every
file to local disk on that guess. A wrong guess costs a full plaintext copy on disk
until the real decision runs; only the later, non-predicting `db.should_qc()` call
is the binding decision.

- [Concepts](#concepts)
- [The pipeline](#the-pipeline)
- [Detailed QC: prediction and decision](#detailed-qc-prediction-and-decision)
- [Parallel instances](#parallel-instances)
- [Failure handling](#failure-handling)
- [Recovery and reruns](#recovery-and-reruns)
- [CLI options](#cli-options)
- [Config keys](#config-keys)

## Concepts

`grzctl process` replaces the individual step-by-step subcommands (`download`,
`decrypt`, `validate`, `encrypt`, `archive`) with one streaming pipeline that avoids
materialising intermediate files on disk for the main pass.

The _interrogation bucket_ is a staging area: each file is uploaded there first,
under its archive key, and only copied to the final archive bucket
once the whole submission has passed. If processing fails, the staged files are
either cleaned up or kept, depending on `archives.interrogation.keep_failed`.

The [README](../README.md#s3-permissions-for-grzctl-process) lists the S3
permissions that `process` needs.

A basic invocation:

```console
$ grzctl --config $CONFIG_PATH process \
    --submission-id $submission_id \
    --output-dir ./out
```

## The pipeline

```
1. download metadata.json from inbox
   add + populate submission's DB row, and record its inbox;
   a case link set by `db case relink` is kept;
   tanG already held by another submission?
   state = ERROR (duplicate_tang); STOP
   case already has a QC-passed initial submission?
   basic_qc_passed = false, state = ERROR (duplicate_initial); STOP
                     │
                     ▼
2. predict detailed QC (detailed_qc.target_percentage > 0 only)
   db.should_qc(predict=True): skips the basic-QC check, stores nothing,
   returns an already-stored decision if there is one
                     │
                     ▼
3. main pass, per file, in a thread pool (--threads files at once),
   skipping outputs that are recorded and still there (see Recovery and reruns):
     inbox object ─► decrypt ─► [predicted yes: write decrypted file to
                                 detailed_qc.local_storage/<submission_id>/files/...]
                  ─► validate (checksum; FASTQ/BAM format)
                  ─► re-encrypt (consented/non-consented key, by consent at the submission date)
                  ─► upload to interrogation bucket (archive key; multipart upload)
   each staged file is recorded in logs/progress_staging.cjson,
   each local copy in logs/progress_local.cjson
                     │
                     ▼
              any file failed? ───────────────────────────────────┐
                     │ no                                          │ yes
                     ▼                                             ▼
5. basic_qc_passed = true                           4. delete prefetched local copies;
   (submission enters the QC queue)                     delete staged objects from
                     │                                   interrogation bucket (unless
                     ▼                                   keep_failed); state = ERROR
6. decide detailed QC                                    STOP
   (target_percentage > 0 only)
   db.should_qc(): reads submitter's QC queue and
   stores selected_for_qc, one decision at a time
   (database-wide lock)
                     │
        ┌────────────┴─────────────┐
        ▼ selected                 ▼ not selected
   QC pass: download/decrypt/     delete prefetched local copies
   checksum-check only files
   without a local copy
   into local storage; write
   local_storage/<submission_id>/
   metadata/metadata.json;
   detailed_qc.auto_run: run
   detailed_qc.shell_command
        └────────────┬─────────────┘
                     ▼
7. upload redacted metadata.json + (redacted) log files to interrogation bucket
                     ▼
8. copy all staged objects to target archive bucket, delete from interrogation bucket
                     ▼
        state = PROCESSED
                     ▼
9. generate Prüfbericht, optionally save (--save-pruefbericht) and
   submit (unless --no-submit-pruefbericht or already REPORTED, with retries)
   (DB states REPORTING → REPORTED; ERROR if BfArM never accepts it, STOP)
                     ▼
10. --clean-inbox (default on): remove submission from inbox
    (DB states CLEANING → CLEANED)
```

Steps 1 to 8 run inside the DB state transition `PROCESSING → PROCESSED` (or `ERROR` on
failure), except the download of `metadata.json`. Processing can create the submission in
the database, so a mistyped submission ID has to fail before that. Steps 9 and 10 each
record their own states, so a failure in them leaves `PROCESSED` in place. A Prüfbericht
that BfArM does not accept keeps the submission in the inbox.

## Detailed QC: prediction and decision

Detailed QC only runs when `detailed_qc.target_percentage` is greater than `0`.

**Prediction (step 2)** happens before the main pass, so the pipeline can decide
whether to prefetch a decrypted copy of each file while it is already streaming
through decrypt anyway. It calls `db.should_qc(predict=True)`, which:

- skips the `basic_qc_passed` check (the submission has not passed it yet),
- stores nothing,
- returns an already-stored `selected_for_qc` decision if one already exists.

The prediction cannot foresee a random selection, because the submission is not in
the QC queue yet (it only enters the queue once `basic_qc_passed` is set to `true`
in step 5). It is a heuristic that is right most of the time. If the DB lookup
fails, the prediction is "no" and a warning is logged.

**Decision (step 6)** happens after `basic_qc_passed = true` is stored, so the
submission is now in the submitter's QC queue. `db.should_qc()` reads that queue and
stores `selected_for_qc` in one transaction, under a database-wide lock:

- PostgreSQL: an advisory transaction lock.
- SQLite: `BEGIN IMMEDIATE`.

So decisions are serialized, also across processes.

- **Selected**: the QC pass downloads, decrypts, and checksum-checks into local storage
  only the files without a local copy (none, if the prediction was right). It then writes
  `<local_storage>/<submission_id>/metadata/metadata.json`. If `detailed_qc.auto_run`
  is true, it runs `detailed_qc.shell_command`. A command that exits non-zero fails the
  run, with the failure reason `detailed_qc_error`.
- **Not selected**: the prefetched local copies are deleted.
- If a later step fails after a positive decision, the local QC data is kept; a
  rerun does not write it again.

## Parallel instances

Two `grzctl process` runs, for different submissions of the same submitter in the
same month, where nobody has been selected for QC yet:

```
 instance A                       DB                        instance B
     │                             │                             │
     │──predict should_qc()──────►│                             │
     │◄──── yes (heuristic) ──────│                             │
     │  prefetch files (assuming the prediction holds)           │
     │                             │◄─────predict should_qc()───│
     │                             │────── yes (heuristic) ────►│
     │                             │     prefetch files (assuming the prediction holds)
     │                             │                             │
     │──basic_qc_passed = true───►│                             │
     │──should_qc(): acquire lock►│                             │
     │   read queue: nobody selected yet                        │
     │   store selected_for_qc = true                            │
     │◄──── true, release lock ───│                             │
     │  selected: keep local copy, run detailed QC               │
     │                             │◄────basic_qc_passed = true─│
     │                             │◄────should_qc(): wait for lock
     │                             │     (waits until A releases it)
     │                             │     read queue: A already selected
     │                             │     store selected_for_qc = false
     │                             │───── false, release lock ─►│
     │                             │   not selected: delete local copies
```

Both instances predict "yes" and both prefetch. The two decisions still run one
after the other, because of the lock: A stores `selected_for_qc = true` first. B's
decision then sees that selection, so the month rule no longer selects B; B is only
selected if another rule still applies (the quarter ratio is at or below target, or B
is the random pick of its block). If B is not selected, it deletes its local copies.

The cost of a wrong "yes" guess is a full plaintext copy on local disk until the
decision runs. The lock is database-wide, not per submitter, because a decision
itself takes only milliseconds.

## Failure handling

If any file fails in the main pass (step 3), the pipeline:

1. deletes the prefetched local copies, if a prediction had written them,
2. deletes the staged objects from the interrogation bucket, unless
   `archives.interrogation.keep_failed` is `true`,
3. fails the run; the DB state becomes `ERROR`.

Step 5 can fail the same way. If another initial submission of
the same case passed basic QC during the main pass, the database rejects
`basic_qc_passed = true`. The pipeline then stores `basic_qc_passed = false` and
handles the failure as above, with the failure reason `duplicate_initial`.

[Error handling](error-handling.md) lists every failure reason and who acts on it.

## Recovery and reruns

Two progress logs live under `<output-dir>/logs/`:

| Log file                 | Written during                                               | What a rerun does with it                                                                                                                                                                          |
| ------------------------ | ------------------------------------------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `progress_staging.cjson` | Main pass (step 3)                                           | Skips validating and staging a file whose re-encrypted copy is recorded and still in the interrogation bucket. The entry also keeps the file's read counts for the read-pair check of its partner. |
| `progress_local.cjson`   | Main pass (with a "yes" prediction) and the QC pass (step 6) | Skips writing a file whose decrypted copy is recorded and still on local storage.                                                                                                                  |

Re-running the command after a partial failure is therefore idempotent at the file
level: each pass downloads a file only for the outputs it is missing, and leaves the
outputs that are still there untouched. A file whose staged copy is gone, for example
after a failed run with `keep_failed: false`, is validated and staged again. A file
whose local copy is gone is written again, by the main pass if the prediction is
"yes", otherwise by the QC pass.

Each interrupt stops a run one level further:

1. The first Ctrl-C starts no queued file and lets the streaming files finish.
2. A second Ctrl-C stops the streaming files at their next chunk.
   Their uploads to the interrogation bucket are aborted, and the progress logs record them as failed.
3. A third Ctrl-C stops waiting for them. The process still exits only once they have stopped.

SIGTERM starts at the second level, because a supervisor follows it with SIGKILL after a grace period.
Either way, the run records `interrupted`.
It leaves the staged and local copies of the finished files in place, so the rerun reuses them.

If the target archive already holds the submission, a rerun streams and copies nothing.
This is the case after a failed Prüfbericht, for example.
The run checks the target archive for the submission's `metadata.json`, which step 8 copies last.
If it is there, the run logs a warning and records `PROCESSED`.
If the database already records `REPORTED`, step 9 does not submit the Prüfbericht again, but still generates and saves it.

Once a cleaning has started, the inbox holds a `cleaning` or `cleaned` marker.
A rerun of `process` then fails with `submission_cleaned`.
`grzctl clean` finishes the cleaning instead, also one that stopped after its deletes.

## CLI options

| Option                                               | Default                 | Effect on this flow                                                                                                                                                                                                                                                  |
| ---------------------------------------------------- | ----------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--threads`                                          | `min(cpu_count, 4)`     | Number of files processed concurrently in the thread pool (step 3 and the QC pass).                                                                                                                                                                                  |
| `--concurrent-uploads`                               | `4`                     | Maximum concurrent part uploads per file's multipart upload to the interrogation bucket. See [Memory use](#memory-use).                                                                                                                                              |
| `--inbox`                                            | `None`                  | Selects which inbox to read from. Defaults to the inbox recorded in the database, the submitter's only inbox, or the only inbox that holds the submission.                                                                                                           |
| `--clean-inbox` / `--no-clean-inbox`                 | `--clean-inbox`         | Whether step 10 removes the submission from the inbox after success.                                                                                                                                                                                                 |
| `--submit-pruefbericht` / `--no-submit-pruefbericht` | `--submit-pruefbericht` | Submits the generated Prüfbericht to BfArM after processing. A server error, a timeout or no connection is retried with backoff.                                                                                                                                     |
| `--save-pruefbericht PATH`                           | `None`                  | Also writes the generated Prüfbericht to `PATH`. A copy with redacted TAN always goes to `logs/pruefbericht.json`.                                                                                                                                                   |
| `--redact-pruefbericht` / `--no-redact-pruefbericht` | `--redact-pruefbericht` | Whether the TAN is redacted in the file written by `--save-pruefbericht`.                                                                                                                                                                                            |

### Memory use

The part buffers of the upload to the interrogation bucket take most of the memory. For one file, they hold up to

```
(C + 1) × P
```

- C is `--concurrent-uploads`.
- P is the part size: `archives.interrogation.s3.multipart_chunksize`, raised to the file size / 1000 for a file
  above 1000 parts.

Decryption, validation and re-encryption add about 60 MiB per file, measured for a gzipped FASTQ. The validators
queue Crypt4GH segments of 64 KiB, and grz-check reads through buffers of 8 MiB.

A run holds this for `--threads` files at once. With the defaults (C = 4, P = 256 MiB, `--threads 4`), that is
about 1.3 GiB per file and 5 GiB per run. A larger file raises P: a 500 GiB file holds 5 parts of 512 MiB.

A future `--max-memory` option should derive `--threads` and the part size from a memory limit.

## Config keys

| Key                                                                                          | Default                               | Effect on this flow                                                                                                                                                          |
| -------------------------------------------------------------------------------------------- | ------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `detailed_qc.local_storage`                                                                  | required, no default                  | Base path for prefetched files and, if selected, the QC pass: `<local_storage>/<submission_id>/files/...` and `/metadata/metadata.json`.                                     |
| `detailed_qc.salt`                                                                           | required, no default                  | Salt for the deterministic random-selection part of `db.should_qc()`.                                                                                                        |
| `detailed_qc.target_percentage`                                                              | `2.0`                                 | Target percentage of a submitter's submissions selected for detailed QC. `0` disables prediction and decision entirely.                                                      |
| `detailed_qc.auto_run`                                                                       | `false`                               | If `true`, runs `detailed_qc.shell_command` right after a selected submission's QC pass completes.                                                                           |
| `detailed_qc.shell_command`                                                                  | a `nextflow run main.nf ...` template | Shell command run when `auto_run` is `true`; templated with `{submission_basepath}`, `{output_basepath}`, `{submission_id}`.                                                 |
| `archives.interrogation`                                                                     | required, no default                  | S3 connection and bucket of the staging area. Every grzctl command loads the whole config, so a config without this section fails for every command, not only for `process`. |
| `archives.interrogation.keep_failed`                                                         | `false`                               | If `true`, leaves a failed submission's staged files in the interrogation bucket instead of deleting them.                                                                   |
| `archives.interrogation.s3.multipart_chunksize`                                              | 256 MiB                               | Preferred part size for uploads to the interrogation bucket. Raised automatically when a file would need more than 1000 parts.                                               |
| `archives.consented.s3.multipart_chunksize`, `archives.non_consented.s3.multipart_chunksize` | 256 MiB                               | Preferred part size for the copy into the final archive. Raised automatically when a file would need more than 1000 parts.                                                   |
