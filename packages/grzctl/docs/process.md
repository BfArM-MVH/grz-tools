# Processing a submission

`grzctl process` runs a submission through the whole pipeline in a single streaming
pass: download the metadata, validate and re-encrypt each file, stage the result in
the *interrogation bucket*, decide whether the submission goes through detailed QC,
then commit everything to the target archive bucket and clean up the inbox. With
`--update-db` it also drives the submission's database state and its QC selection.

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

The *interrogation bucket* is a staging area: each file is uploaded there first,
under its archive key, and only copied to the final archive bucket
once the whole submission has passed. If processing fails, the staged files are
either cleaned up or kept, depending on `archives.interrogation.keep_failed`.

A basic invocation:

```console
$ grzctl --config $CONFIG_PATH process \
    --submission-id $submission_id \
    --output-dir ./out \
    --update-db
```

## The pipeline

```
1. download metadata.json from inbox
   --update-db: add + populate submission's DB row;
                case already has a QC-passed initial submission?
                basic_qc_passed = false, state = ERROR (duplicate_initial); STOP
                     │
                     ▼
2. predict detailed QC (--update-db and detailed_qc.target_percentage > 0 only)
   db.should_qc(predict=True): skips the basic-QC check, stores nothing,
   returns an already-stored decision if there is one
                     │
                     ▼
3. main pass, per file, in a thread pool (--threads files at once):
     inbox object ─► decrypt ─► [predicted yes: write decrypted file to
                                 detailed_qc.local_storage/<submission_id>/files/...]
                  ─► validate (checksum; FASTQ/BAM format)
                  ─► re-encrypt (consented/non-consented key, by consent at run time)
                  ─► upload to interrogation bucket (archive key; multipart upload)
   each file that succeeds is recorded in logs/progress_processing.cjson;
   with the prediction, also in logs/progress_qc.cjson
                     │
                     ▼
              any file failed? ───────────────────────────────────┐
                     │ no                                          │ yes
                     ▼                                             ▼
5. --update-db: basic_qc_passed = true              4. delete prefetched local copies +
   (submission enters the QC queue)                     progress_qc.cjson; delete staged
                     │                                   objects from interrogation bucket
                     ▼                                   (unless keep_failed); state = ERROR
6. decide detailed QC (--update-db and                   with --update-db; STOP
   target_percentage > 0 only)
   db.should_qc(): reads submitter's QC queue and
   stores selected_for_qc, one decision at a time
   (database-wide lock)
                     │
        ┌────────────┴─────────────┐
        ▼ selected                 ▼ not selected
   QC pass: download/decrypt/     delete prefetched local copies
   checksum-check only files      + progress_qc.cjson
   not yet in progress_qc.cjson
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
9. --clean-inbox (default on): remove submission from inbox
   (DB states CLEANING → CLEANED)
                     ▼
        state = PROCESSED (--update-db)
                     ▼
   generate Prüfbericht, optionally save (--save-pruefbericht) and
   submit (--submit-pruefbericht, with retries)
```

The duplicate-initial check in step 1 and steps 2-9 run inside the DB state transition `PROCESSING → PROCESSED` (or `ERROR`
on failure) when `--update-db` is set.

## Detailed QC: prediction and decision

Detailed QC only runs when `--update-db` is set and `detailed_qc.target_percentage`
is greater than `0`. With `--no-update-db` there is no detailed QC at all: the
selection itself writes to the database, so it is skipped outright.

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

- **Selected**: the QC pass downloads, decrypts, and checksum-checks only the files
  not yet recorded in `progress_qc.cjson` (none, if the prediction was right and
  nothing was skipped) into local storage. It then writes
  `<local_storage>/<submission_id>/metadata/metadata.json`. If `detailed_qc.auto_run`
  is true, it runs `detailed_qc.shell_command`.
- **Not selected**: the prefetched local copies and `progress_qc.cjson` are deleted.
- If a later step fails after a positive decision, the local QC data is kept; a
  rerun skips the files already recorded in `progress_qc.cjson`.

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
     │                             │                + progress_qc.cjson
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

1. deletes the prefetched local copies and `progress_qc.cjson`, if a prediction had
   written them,
2. deletes the staged objects from the interrogation bucket, unless
   `archives.interrogation.keep_failed` is `true`,
3. fails the run; with `--update-db`, the DB state becomes `ERROR`.

With `--update-db`, step 5 can fail the same way. If another initial submission of
the same case passed basic QC during the main pass, the database rejects
`basic_qc_passed = true`. The pipeline then stores `basic_qc_passed = false` and
handles the failure as above, with the failure reason `duplicate_initial`.

## Recovery and reruns

Three progress logs live under `<output-dir>/logs/`:

| Log file | Written during | What a rerun does with it |
| --- | --- | --- |
| `progress_download.cjson` | Metadata download | Tracks the metadata download itself. |
| `progress_processing.cjson` | Main pass (step 3) | Files already recorded as processed are skipped, and are not written to local storage again. |
| `progress_qc.cjson` | Main pass (with a "yes" prediction) and the QC pass (step 6) | The QC pass only re-downloads and re-checks files not yet recorded here. |

Re-running the command after a partial failure is therefore idempotent at the file
level: already-completed files are skipped in both the main pass and, if selected,
the QC pass. `progress_qc.cjson` is deleted whenever the submission's local QC copy
is deleted (failure, or a "not selected" decision), so a later rerun starts QC
prefetching from scratch if it is picked again.

## CLI options

| Option | Default | Effect on this flow |
| --- | --- | --- |
| `--update-db` / `--no-update-db` | `--update-db` | Turns the DB row, `basic_qc_passed`, the detailed-QC prediction and decision, and the `PROCESSING`/`PROCESSED`/`ERROR` state transitions on or off. With `--no-update-db`, detailed QC never runs, and duplicate initial submissions are not detected. |
| `--threads` | `min(cpu_count, 4)` | Number of files processed concurrently in the thread pool (step 3 and the QC pass). |
| `--concurrent-uploads` | `4` | Maximum concurrent part uploads per file's multipart upload to the interrogation bucket. |
| `--inbox-bucket` | `None` | Selects which inbox to read from, if the submitter has more than one configured. |
| `--clean-inbox` / `--no-clean-inbox` | `--clean-inbox` | Whether step 9 removes the submission from the inbox after success. |
| `--submit-pruefbericht` / `--no-submit-pruefbericht` | `--no-submit-pruefbericht` | Submits the generated Prüfbericht to BfArM after processing, with retries. |
| `--save-pruefbericht PATH` | `None` | Also writes the generated Prüfbericht to `PATH`. A copy with redacted TAN always goes to `logs/pruefbericht.json`. |
| `--redact-pruefbericht` / `--no-redact-pruefbericht` | `--redact-pruefbericht` | Whether the TAN is redacted in the file written by `--save-pruefbericht`. |

## Config keys

| Key | Default | Effect on this flow |
| --- | --- | --- |
| `detailed_qc.local_storage` | required, no default | Base path for prefetched files and, if selected, the QC pass: `<local_storage>/<submission_id>/files/...` and `/metadata/metadata.json`. |
| `detailed_qc.salt` | required, no default | Salt for the deterministic random-selection part of `db.should_qc()`. |
| `detailed_qc.target_percentage` | `2.0` | Target percentage of a submitter's submissions selected for detailed QC. `0` disables prediction and decision entirely. |
| `detailed_qc.auto_run` | `false` | If `true`, runs `detailed_qc.shell_command` right after a selected submission's QC pass completes. |
| `detailed_qc.shell_command` | a `nextflow run main.nf ...` template | Shell command run when `auto_run` is `true`; templated with `{submission_basepath}`, `{output_basepath}`, `{submission_id}`. |
| `archives.interrogation.keep_failed` | `false` | If `true`, leaves a failed submission's staged files in the interrogation bucket instead of deleting them. |
| `archives.interrogation.s3.multipart_chunksize` | 256 MiB | Preferred part size for uploads to the interrogation bucket. Raised automatically when a file would need more than 1000 parts. |
