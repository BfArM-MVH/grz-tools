# Error handling

When a grzctl step fails, the database records the submission state `ERROR` together
with a _failure reason_. The reason names the cause of the failure, and with it who
has to act: the submitter, the operator, or the GRZ itself. Where the failure happened
goes into the state's data.

- [What a failure records](#what-a-failure-records)
- [Failure reasons](#failure-reasons)
- [Several failures in one run](#several-failures-in-one-run)
- [Raising errors](#raising-errors)
- [S3 errors](#s3-errors)

## What a failure records

Each step runs in a database context. The context records the step's start state, such
as `PROCESSING`, when the step begins. It records the end state when the step succeeds.
When the step fails, it records `ERROR` with:

- `failure_reason`: one of the [failure reasons](#failure-reasons).
- `data.error`: the message of the failure that decided the reason.
- `data.errors`: one entry per failed file, with `file`, `reason` and `message`. Only a
  step that streams files writes it.

Each step records its own outcome, and records it once. `grzctl process` runs three
steps, one after another:

| Step       | States                     | Covers                                                                                                          |
| ---------- | -------------------------- | --------------------------------------------------------------------------------------------------------------- |
| Processing | `PROCESSING` → `PROCESSED` | metadata download and parsing, database population, duplicate checks, file streaming, detailed QC, archive copy |
| Reporting  | `REPORTING` → `REPORTED`   | generating and submitting the Prüfbericht                                                                       |
| Cleaning   | `CLEANING` → `CLEANED`     | removing the submission from the inbox                                                                          |

A failed step ends the run, so the inbox keeps the submission until BfArM has accepted
its Prüfbericht.

A failure before a step's context is open reaches only the console and the log. Examples
are an unreachable database and a missing BfArM credential.

## Failure reasons

| Reason                  | Who acts  | Meaning                                                                                                                                                                        |
| ----------------------- | --------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `file_not_found`        | submitter | The inbox lacks the metadata or a file that the metadata lists.                                                                                                                |
| `validation_error`      | submitter | The metadata or the content of a file breaks the specification: schema, checksum, file format, or the read counts of a read pair.                                              |
| `decryption_error`      | submitter | A file cannot be decrypted with the GRZ private key. If every submission fails this way, check the configured private key.                                                     |
| `duplicate_tang`        | submitter | Another submission already used the tanG.                                                                                                                                      |
| `duplicate_initial`     | submitter | The case already has an initial submission that passed basic QC.                                                                                                               |
| `incomplete_submission` | operator  | A step-by-step command ran before an earlier step had passed for every file.                                                                                                   |
| `interrupted`           | operator  | Ctrl-C or SIGTERM stopped the run. A rerun resumes it.                                                                                                                         |
| `configuration_error`   | GRZ       | The GRZ setup is wrong or incomplete: a key or credential is missing, S3 rejects the credentials, a configured bucket does not exist, or BfArM refuses the client credentials. |
| `transfer_error`        | GRZ       | Moving data to or from S3 failed, or S3 stored other bytes than were sent. A rerun usually succeeds. If the failure repeats, the message names the S3 error code.              |
| `encryption_error`      | GRZ       | Re-encrypting a file for the archive failed.                                                                                                                                   |
| `detailed_qc_error`     | GRZ       | The detailed QC workflow exited with an error.                                                                                                                                 |
| `reporting_error`       | GRZ       | The Prüfbericht could not be generated, or BfArM did not accept it. `grzctl process` retries a submission that fails this way.                                                 |
| `unknown`               | GRZ       | No code path expected this error. Treat it as a bug.                                                                                                                           |

Older states can carry `network_error` or `upload_error`. Both mean what
`transfer_error` means. grzctl no longer writes them, and `grzctl db submission update`
does not offer them.

## Several failures in one run

`grzctl process` streams files in parallel, so one run can collect several file errors.
The most decisive of them sets the reason, in this order:

1. a submitter reason
2. `configuration_error`
3. `transfer_error`
4. any other GRZ reason

A rejected submission stays rejected after any rerun, so its reason decides what happens
next. Within the same rank, the first failed file in metadata order wins. `data.errors`
lists every file error, whatever its reason.

## Raising errors

Every expected failure is a `GrzError`. Its class sets the failure reason:

```
GrzError
├── SubmissionRejectedError
│   ├── MissingSubmissionFileError   file_not_found
│   ├── SubmissionValidationError    validation_error
│   ├── DecryptionError              decryption_error
│   └── DuplicateUploadError         duplicate_tang
├── IncompleteSubmissionError        incomplete_submission
├── ConfigurationError               configuration_error
├── TransferError                    transfer_error
│   ├── DownloadError
│   │   └── MissingObjectError
│   └── UploadError
│       └── UploadIntegrityError
├── EncryptionError                  encryption_error
├── DetailedQCError                  detailed_qc_error
└── ReportingError                   reporting_error
```

The errors of grz-db that reject a duplicate map to `duplicate_tang` and
`duplicate_initial`. A `KeyboardInterrupt` maps to `interrupted`, and grzctl turns
SIGTERM into one. Every other exception maps to `unknown`.

When you raise an error:

- Wrap a library error once, at the code that knows what the call meant. An S3 client
  does not know that a bucket is the inbox, but the processor does. So the processor
  turns a missing inbox object into a `MissingSubmissionFileError`.
- Chain the cause with `raise ... from e`. The log then keeps the library's traceback.
- Do not raise a builtin exception, such as `ValueError`, for an expected failure. Its
  reason is `unknown`.
- Do not call `sys.exit` below the command line layer. `SystemExit` passes every
  `except Exception`, so no cleanup runs.
- Give every new `GrzError` subclass a failure reason. A test fails for a subclass that
  maps to `unknown`.
- Cleanup after a failure, such as deleting staged objects, logs its own errors. It never
  replaces the original error.

## S3 errors

The S3 boundary sorts the error codes like this:

| S3 answer                                                     | Raised as                    | Reason                |
| ------------------------------------------------------------- | ---------------------------- | --------------------- |
| 404 for an inbox object                                       | `MissingSubmissionFileError` | `file_not_found`      |
| 404 for any other object                                      | `MissingObjectError`         | `transfer_error`      |
| `InvalidAccessKeyId`, `SignatureDoesNotMatch`, `NoSuchBucket` | `ConfigurationError`         | `configuration_error` |
| `AccessDenied` and anything else                              | `TransferError`              | `transfer_error`      |

Missing credentials count as a configuration error as well.

`AccessDenied` is not a configuration error. S3 also answers 403 for a missing object if
the credentials lack the permission to list the bucket.
