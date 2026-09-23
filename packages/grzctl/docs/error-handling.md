# Error handling

When a grzctl step fails, the database records the submission state `ERROR` together with a _failure reason_.
The reason names the cause of the failure, and with it who has to act.
That is either the LE (Leistungserbringer), which sent the submission, or the GRZ, which runs grzctl.
Where the failure happened goes into the state's data.

- [What a failure records](#what-a-failure-records)
- [Failure reasons](#failure-reasons)
- [Raising errors](#raising-errors)
- [S3 errors](#s3-errors)
- [BfArM errors](#bfarm-errors)

## What a failure records

Each grzctl command that takes `--update-db` runs its step in a database context.
The context records the step's start state, such as `DOWNLOADING`, when the step begins.
It records the end state when the step succeeds.
When the step fails, it records `ERROR` with:

- `failure_reason`: one of the [failure reasons](#failure-reasons).
- `data.error`: the message of the failure that decided the reason.

Each step records its own outcome, and records it once.

A failure before a step's context is open reaches only the console and the log.
An example is an unreachable database.
A step creates a missing submission in the database only if it starts at `PROCESSING` or `UPLOADING`.
Every other step stops for a mistyped submission ID before its context is open, so the ID leaves no trace in the database.

## Failure reasons

| Reason                          | Who acts | Meaning                                                                                                                                       |
| ------------------------------- | -------- | --------------------------------------------------------------------------------------------------------------------------------------------- |
| `file_not_found`                | LE       | The inbox lacks the metadata or a file that the metadata lists.                                                                               |
| `validation_error`              | LE       | The metadata or the content of a file breaks the specification: schema, checksum, file format, or the read counts of a read pair.             |
| `decryption_error`              | LE       | A file cannot be decrypted with the GRZ private key. If every submission fails this way, check the configured private key.                    |
| `duplicate_tang`                | LE       | Another submission already used the tanG.                                                                                                     |
| `duplicate_initial`             | LE       | The case already has an initial submission that passed basic QC.                                                                              |
| `incomplete_submission`         | GRZ      | A command ran before an earlier step had passed for every file.                                                                               |
| `submission_cleaned`            | GRZ      | `grzctl clean` has started on the submission, so the inbox no longer holds it.                                                                |
| `interrupted`                   | GRZ      | Ctrl-C or SIGTERM stopped the run. A rerun resumes it.                                                                                        |
| `configuration_error`           | GRZ      | The setup is wrong or incomplete: a key or credential is missing, S3 or BfArM rejects the credentials, or a configured bucket does not exist. |
| `transfer_error`                | GRZ      | Moving data to or from S3 or BfArM failed. A rerun usually succeeds. If the failure repeats, the message names the error code or HTTP status. |
| `encryption_error`              | GRZ      | Re-encrypting a file for the archive failed.                                                                                                  |
| `detailed_qc_error`             | GRZ      | The detailed QC workflow exited with an error.                                                                                                |
| `pruefbericht_generation_error` | GRZ      | The Prüfbericht cannot be generated, because the database lacks the submission or data that the Prüfbericht needs.                            |
| `pruefbericht_rejected`         | GRZ      | BfArM rejected the Prüfbericht. Sending the same one again fails again.                                                                       |
| `unknown`                       | GRZ      | No code path accounts for this error. Treat it as a bug.                                                                                      |

No grzctl step raises the error behind `detailed_qc_error` yet.
`grzctl db submission update --failure-reason` can record it.

Older states can carry `network_error` or `upload_error`.
Both mean what `transfer_error` means.
grzctl no longer writes them, and `grzctl db submission update` does not offer them.

## Raising errors

Every expected failure is a `GrzError`.
Its class sets the failure reason:

```
GrzError
├── SubmissionRejectedError
│   ├── MissingSubmissionFileError   file_not_found
│   ├── SubmissionValidationError    validation_error
│   ├── DecryptionError              decryption_error
│   └── DuplicateUploadError         duplicate_tang
├── IncompleteSubmissionError        incomplete_submission
├── SubmissionCleanedError           submission_cleaned
├── ConfigurationError               configuration_error
├── TransferError                    transfer_error
│   ├── DownloadError
│   │   └── MissingObjectError
│   ├── UploadError
│   └── NetworkError
├── EncryptionError                  encryption_error
├── DetailedQCError                  detailed_qc_error
├── PruefberichtGenerationError      pruefbericht_generation_error
└── PruefberichtRejectedError        pruefbericht_rejected
```

The errors of grz-db that reject a duplicate map to `duplicate_tang` and `duplicate_initial`.
A `KeyboardInterrupt` maps to `interrupted`, and grzctl turns SIGTERM into one.
Every other exception maps to `unknown`.

When you raise an error:

- Wrap a library error once, at the code that knows what the call meant.
  An S3 client does not know that a bucket is the inbox, but the download worker does.
  So the download worker turns a missing inbox object into a `MissingSubmissionFileError`.
- Chain the cause with `raise ... from e`.
  The log then keeps the library's traceback.
- Do not raise a builtin exception, such as `ValueError`, for an expected failure.
  Its reason is `unknown`.
- Do not call `sys.exit` below the command line layer.
  `SystemExit` passes every `except Exception`, so no cleanup runs.
- Give every new `GrzError` subclass a failure reason.
  A test fails for a subclass that maps to `unknown`.
- Cleanup after a failure logs its own errors.
  It never replaces the original error.

## S3 errors

The S3 boundary sorts the error codes like this:

| S3 answer                                                     | Raised as                    | Reason                |
| ------------------------------------------------------------- | ---------------------------- | --------------------- |
| `NoSuchKey` for an inbox object                               | `MissingSubmissionFileError` | `file_not_found`      |
| `NoSuchKey` for any other object                              | `MissingObjectError`         | `transfer_error`      |
| `InvalidAccessKeyId`, `SignatureDoesNotMatch`, `NoSuchBucket` | `ConfigurationError`         | `configuration_error` |
| `AccessDenied` and anything else                              | `TransferError`              | `transfer_error`      |

Missing credentials count as a configuration error as well, and so does a bucket name that botocore rejects before it sends a request.
The upload worker takes the error code from the `ClientError` that `S3Transfer` wraps in `S3UploadFailedError`.

S3 answers a HEAD request without a body, so botocore reports only the HTTP status as the error code.
A `403` can then mean rejected credentials, and a `404` a missing bucket.
For these two codes, `head_object()` sends a GET request for the first byte of the object and sorts the error code of that answer.
Downloads and uploads start with `head_object()`, so their first request already tells a faulty setup from a missing object.

`grzctl clean` marks its start with a `<submission_id>/cleaning` object and its end with `<submission_id>/cleaned`.
It also empties the metadata.json, whose `LastModified` then is the time of cleaning.
So before grzctl downloads the metadata.json or reads its upload date, it checks for both markers and for an empty metadata.json.
If it finds one of them, it raises a `SubmissionCleanedError`.

`AccessDenied` is not a configuration error.
S3 also answers 403 for a missing object if the credentials lack the permission to list the bucket.
For the same reason, `AccessDenied` for the submission's metadata does not tell whether the bucket holds the submission.
An upload or an archival then logs a warning that it cannot check for an earlier run, and it continues.

## BfArM errors

`grzctl pruefbericht submit` sorts the answers of BfArM like this:

| BfArM answer                                 | Raised as                   | Reason                  |
| -------------------------------------------- | --------------------------- | ----------------------- |
| `401` or `403`                               | `ConfigurationError`        | `configuration_error`   |
| Any other client error for the token request | `ConfigurationError`        | `configuration_error`   |
| Any other client error for the Prüfbericht   | `PruefberichtRejectedError` | `pruefbericht_rejected` |
| A server error, a timeout, or no connection  | `NetworkError`              | `transfer_error`        |

The token request carries only the client credentials, so any client error means that they are wrong.
A `401` for the Prüfbericht first makes grzctl request a new token and send the Prüfbericht once more.
