# Registering change requests

A *change request* records that a submission has to be modified, deleted, or
transferred — together with an audit trail of **who** asked for it, **when**,
and **the content** of their request. This page walks through the full CLI
workflow.

> Change requests are rare but consequential: each one is signed by the acting
> author and stored permanently. The CLI therefore enforces that the audit
> fields are present before anything is written.

- [Concepts](#concepts)
- [The three change types](#the-three-change-types)
- [Quick start](#quick-start)
- [Getting a template](#getting-a-template)
- [Supplying the data](#supplying-the-data)
  - [`--data-file` (YAML or JSON file)](#--data-file-yaml-or-json-file)
  - [`--data` (inline JSON)](#--data-inline-json)
- [Attaching the original request (`--raw-content`)](#attaching-the-original-request---raw-content)
- [Required fields and validation](#required-fields-and-validation)
  - [`change-request-validate` (offline, no config)](#change-request-validate--offline-no-config)
  - [`--dry-run` (online, needs config)](#--dry-run--online-needs-config)
- [Extra / type-specific data](#extra--type-specific-data)
- [Inspecting stored change requests](#inspecting-stored-change-requests)
- [Runnable end-to-end demo](#runnable-end-to-end-demo)
- [Command reference](#command-reference)

## Concepts

A change request is attached to an existing submission and captures:

| Field | Required | Description |
| --- | --- | --- |
| `requester_name` | yes | Full name of the person who requested the change. |
| `requester_email` | yes | Email address of the requester. |
| `requested_at` | yes | Date the change was requested, as `YYYY-MM-DD`. |
| `request_email_content` | yes\* | Verbatim text of the request (multi-line allowed). |
| `request_raw_content` | no | Optional binary attachment (PDF/PNG), supplied via `--raw-content`. |
| `data` | no | Free-form key/value extras (type-specific notes, ticket IDs, …). |

\* You must provide **either** `request_email_content` **or** `--raw-content`
(you may provide both). This lets you record a request that arrived only as a
signed PDF without also transcribing it.

The audit columns are nullable in the database so that historical rows created
before these columns existed can still be loaded. **Required-ness is enforced by
the CLI**, not the schema — new entries always get the full audit trail.

## The three change types

`CHANGE` is one of (case-insensitive):

- **`Modify`** — the submission's data needs to be corrected.
- **`Delete`** — the submission must be deleted.
- **`Transfer`** — the submission is being transferred elsewhere.

The audit fields are identical across all three — who/what/when is universal.
Only the human-facing template guidance differs per type.

## Quick start

1. Get a template for the change type you need
`grzctl change-request-template Delete > delete-request.yaml`

2. Modify the generated template file with a text editor to fill in needed information

3. Validate the file offline — no config, no DB, cannot write anything
```
grzctl change-request-validate Delete --data-file delete-request.yaml
```

4. Write it for real
```
grzctl db --config-file $CONFIG_PATH submission change-request \
    $submission_id Delete \
    --data-file delete-request.yaml
```

The submission must already exist in the DB. If it does not, add it first:

`grzctl db --config-file $CONFIG_PATH submission add $submission_id`.

## Getting a template

`change-request-template` prints a fill-in YAML skeleton for a change type. It
is a **top-level** command — it needs no `--config-file` or DB credentials,
because it only prints text:

```console
$ grzctl change-request-template Delete
```

```yaml
# Full name of the person requesting the change
requester_name: <FILL IN requester name>

# Email address of the requester
requester_email: <FILL IN email address>

# Date the change was requested (YYYY-MM-DD)
requested_at: YYYY-MM-DD

# Verbatim email text from the requester (multi-line allowed).
# Required unless --raw-content (e.g. PDF) is supplied via the CLI.
request_email_content: |
  <FILL IN paste verbatim email content here — describe what to delete and why>

# Optional type-specific extras (free-form key/value pairs)
data: {}
```

The `<FILL IN …>` placeholders are not just hints: the validator **rejects** any
value still containing the `<FILL IN` marker, so an unedited template cannot be
submitted by accident.

## Supplying the data

There are two mutually exclusive ways to pass the change-request fields:
`--data-file` and `--data`. Passing both is an error.

### `--data-file` (YAML or JSON file)

The recommended path for anything non-trivial, since request emails are usually
multi-line. The file is parsed as YAML, which is a superset of JSON — so `.yaml`,
`.yml`, and `.json` files all work, and the extension is only a hint (a
correctly-formed file loads even if it carries the "wrong" extension). The top
level of the file must be a mapping/object of the change-request fields.

```console
$ grzctl db --config-file $CONFIG_PATH submission change-request \
    $submission_id Delete \
    --data-file delete-request.yaml
```

A filled-in `delete-request.yaml`:

```yaml
requester_name: Erika Mustermann
requester_email: erika.mustermann@example.org
requested_at: 2026-03-27
request_email_content: |
  Liebe Kolleginnen und Kollegen,
  bitte den Datensatz mit der Submission-ID $submission_id löschen.
  Vielen Dank, Erika
data:
  ticket: HELPDESK-1234
```

### `--data` (inline JSON)

Handy for scripting or one-liners where you don't want a file on disk:

```console
$ grzctl db --config-file $CONFIG_PATH submission change-request \
    $submission_id Delete \
    --data '{"requester_name":"Erika Mustermann","requester_email":"erika.mustermann@example.org","requested_at":"2026-03-27","request_email_content":"Bitte löschen."}'
```

The same required-field and validation rules apply as for `--data-file`.

## Attaching the original request (`--raw-content`)

Requests often arrive as a signed PDF or a screenshot. Attach the original file
with `--raw-content`:

```console
$ grzctl db --config-file $CONFIG_PATH submission change-request \
    $submission_id Delete \
    --data-file delete-request.yaml \
    --raw-content signed-deletion-notice.pdf
```

- Supported types: **PDF** and **PNG**.
- The type is inferred from the file extension (`.pdf`, `.png`) and then
  **verified against the file's magic bytes** (`%PDF-` for PDF, the PNG
  signature for PNG). A mismatch or unsupported extension is rejected.
- If you supply `--raw-content`, `request_email_content` becomes optional — the
  attachment can stand in for the transcribed email.

## Required fields and validation

The same input checks apply everywhere a change request is validated:

1. `requester_name`, `requester_email`, and `requested_at` are all present.
2. At least one of `request_email_content` or `--raw-content` is provided.
3. No field still holds a template placeholder — neither the `<FILL IN …>`
   markers nor the `YYYY-MM-DD` date placeholder.
4. `requested_at` is a valid `YYYY-MM-DD` date.
5. The whole record validates against the change-request model (a final backstop
   for anything the checks above don't cover).

If anything is wrong, **all** problems are reported together as a single list —
so you can fix them in one pass rather than one error per run — followed by a
pointer to `grzctl change-request-template <CHANGE>` for a filled-in example. For
example, an unedited template reports each `<FILL IN>` field plus the
`YYYY-MM-DD` date placeholder at once:

```console
$ grzctl change-request-validate Delete --data-file unedited-template.yaml
Error: the change-request input needs attention:
  • requester_name: still has the '<FILL IN ...>' placeholder — replace it
  • requester_email: still has the '<FILL IN ...>' placeholder — replace it
  • request_email_content: still has the '<FILL IN ...>' placeholder — replace it
  • requested_at: still the 'YYYY-MM-DD' placeholder — use a real date, e.g. 2026-03-27
See `grzctl change-request-template Delete` for a filled-in example.
```

There are two ways to run these checks, differing only in whether they touch the
database:

### `change-request-validate` — offline, no config

```console
$ grzctl change-request-validate Delete --data-file delete-request.yaml
```

This is a **top-level** command: it needs no `--config-file`, opens no database
connection, and **cannot write anything**. Use it as the routine "is my file
correct?" check — because there is no DB in the picture, there is no way to
accidentally register a live change request. On success it echoes back the
validated fields (binary content elided as a byte count).

Since change requests are rare and consequential, prefer this offline check
while iterating on your data file, and only reach for the DB command once it
passes.

### `--dry-run` — online, needs config

```console
$ grzctl db --config-file $CONFIG_PATH submission change-request \
    $submission_id Delete \
    --data-file delete-request.yaml --dry-run
```

`--dry-run` runs the identical input checks **and additionally** connects to the
database to confirm the submission exists — but still writes nothing. Use it
when you want that extra "would this actually apply to a real submission?"
confirmation. It requires a valid `--config-file` (including author keys),
because it is a subcommand of the `db` group.

## Extra / type-specific data

Anything beyond the known audit fields is preserved as free-form extras and
stored in the `data` column. You can supply extras either:

- nested under a `data:` key (recommended, explicit), or
- as unknown top-level keys (these are collected into `data` automatically).

```yaml
requester_name: Erika Mustermann
requester_email: erika.mustermann@example.org
requested_at: 2026-03-27
request_email_content: |
  Bitte löschen.
data:
  ticket: HELPDESK-1234
  reason: withdrawn-consent
```

## Inspecting stored change requests

List the change requests recorded for submissions:

```console
$ grzctl db --config-file $CONFIG_PATH list-change-requests
$ grzctl db --config-file $CONFIG_PATH list-change-requests --json
```

The `--json` form includes the audit columns; binary `request_raw_content` is
stored as-is and can be large.

## Runnable end-to-end demo

For a hands-on walkthrough that spins up a throwaway SQLite database and
exercises the template, validation, dry-run, real-write, and `--raw-content`
paths end to end, run:

```console
$ uv run python packages/grzctl/examples/demo_change_request.py
# add --keep to preserve the temp workspace for inspection
```

The demo is kept in lockstep with the CLI (it drives the same Click commands),
so it doubles as living documentation of the flow described above. See
[`examples/demo_change_request.py`](../examples/demo_change_request.py).

## Command reference

| Command | Purpose |
| --- | --- |
| `grzctl change-request-template CHANGE` | Print a fill-in YAML template (no config needed). |
| `grzctl change-request-validate CHANGE` | Validate an input file offline — no config, no DB, no writes. |
| `grzctl db … submission change-request SUBMISSION_ID CHANGE` | Register a change request. |
| `grzctl db … submission add SUBMISSION_ID` | Add a submission (prerequisite). |
| `grzctl db … list-change-requests [--json]` | List stored change requests. |

Options for `submission change-request`:

| Option | Description |
| --- | --- |
| `--data-file PATH` | JSON/YAML file with the change-request fields. Mutually exclusive with `--data`. |
| `--data JSON` | Inline JSON with the change-request fields. Mutually exclusive with `--data-file`. |
| `--raw-content PATH` | Optional PDF/PNG attachment; type verified by magic bytes. |
| `--dry-run` | Validate and check the submission exists, but write nothing. |

`CHANGE` is one of `Modify`, `Delete`, `Transfer` (case-insensitive).
