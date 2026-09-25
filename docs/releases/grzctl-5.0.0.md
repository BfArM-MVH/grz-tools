## Upgrade guide: grzctl 5.0.0

This guide is for GRZ operators. It lists what to do before and after the update, and what behaves differently. The changelog of each package lists every change, see the [releases](https://github.com/BfArM-MVH/grz-tools/releases).

The main changes:

- grzctl reads one config file for all commands (#635).
- grzctl tracks cases with `grzctl db case` (#633), see [Case tracking](../../packages/grzctl/docs/case-tracking.md).
- Detailed QC reports a deviation without failing (#656).
- A failed step records why it failed (#690), see [Error handling](../../packages/grzctl/docs/error-handling.md).
- The grzctl config names each crypt4gh key where grzctl uses it (#694), see [Crypt4GH keys](../../packages/grzctl/docs/crypt4gh-keys.md).
- The database records the inbox of each submission, so `download`, `clean` and `decrypt` find it without `--inbox` (#696).

### Versions

- `grzctl` `5.0.0` (from `4.0.0`)
- `grz-common` `4.0.0` (from `3.0.0`)
- `grz-db` `4.0.0` (from `3.0.0`)
- `grz-pydantic-models` `4.0.0` (from `3.0.0`)
- `grz-pydantic-models-testing` `1.1.0` (from `1.0.0`)
- `grz-cli` `3.0.0` (from `2.0.0`)

---

## Migration instructions

> ⚠️ **Plan for downtime.** Stop all automated worker scripts before you start.

> ⚠️ **Back up the database before `grzctl db upgrade`.** For PostgreSQL, for example: `pg_dump --format=custom --file=grz-db-before-5.0.0.dump <database>`. For SQLite, copy the database file.

### 1. Update to grzctl v5.0.0

Check the versions after the update:

```console
$ grzctl --version
grzctl v5.0.0
grz-common v4.0.0
grz-db v4.0.0
grz-pydantic-models v4.0.0
grz-check v0.4.0
```

grzctl no longer installs grz-cli (#675). If you use grz-cli, for example for `grz-cli upload`, install it as well.

### 2. Write the unified config file (#635)

> ℹ️ This concerns the grzctl config only. The LE configs for grz-cli stay unchanged.

- grzctl reads one file for all commands, by default `~/.config/grzctl/config.yaml`.
- `grzctl --config PATH <command>` replaces `--config-file` on each command. grzctl no longer merges several files.
- Environment variables override the file, for example `GRZ_LEISTUNGSERBRINGER__123456789__INBOX_BUCKETS__MAIN__PRIVATE_KEY_PATH`.

All five top-level sections are required:

```yaml
inbox_defaults: &inbox_defaults # optional; grzctl ignores this section, and the inboxes merge it with <<
  endpoint_url: https://s3.example.org
  private_key_path: /path/to/grz.sec
leistungserbringer:
  "123456789": # LE ID, quoted
    alias: "FOO" # optional
    inbox_buckets:
      main: # the inbox name, which --inbox takes
        <<: *inbox_defaults
        bucket: le-123456789 # optional, defaults to the inbox name
        access_key: ...
        secret: ...
archives:
  consented:
    s3:
      endpoint_url: https://s3.example.org
      bucket: grz-consented
    public_key_path: /path/to/consented.pub
  non_consented:
    s3:
      endpoint_url: https://s3.example.org
      bucket: grz-non-consented
    public_key_path: /path/to/non_consented.pub
db:
  database_url: postgresql+psycopg://...
  author:
    name: ...
    private_key_path: /path/to/author.sec
  known_public_keys: # optional, one "<format> <key> <author name>" entry per key
    - ssh-ed25519 AAAA... alice
  # known_public_keys_file: /path/to/known_public_keys # the same entries as a file, instead of the list
identifiers:
  grz: GRZABC123
pruefbericht:
  authorization_url: https://...
  client_id: ...
  client_secret: ...
  api_base_url: https://...
```

`keys.grz_private_key_path` goes away. Each inbox names its private key, and `encrypt` signs with it (#694, #696). Every key field takes the key inline as `<name>` or a file as `<name>_path`. Setting both is an error. A private key has an optional `<name>_passphrase`. [Crypt4GH keys](../../packages/grzctl/docs/crypt4gh-keys.md) lists every key and its field.

`db.known_public_keys` now lists the keys. Move a path to a file to `db.known_public_keys_file` (#693). If neither is set, grzctl reads `~/.config/grzctl/known_public_keys`.

Check it with `grzctl dump-config`. Add `--reveal-secrets` to show the secrets (#680).

### 3. Upgrade the database

1. Run `grzctl db upgrade`. It adds case tracking (#633), the processing states (#681), the new failure reasons (#690) and the inbox of a submission (#696), and renames `submissions.pseudonym` to `local_case_id`.
2. Run `grzctl db case list-unlinked`, and resolve what it lists as described in [Upgrading an existing database](../../packages/grzctl/docs/case-tracking.md#upgrading-an-existing-database).

### 4. Backfill the lossless metadata (#654)

Backfill the lossless `submission_metadata` into the existing submissions. One run covers both archives:

```bash
grzctl db backfill --dry-run --allow-overwrite submission_metadata   # preview
grzctl db backfill --allow-overwrite submission_metadata
```

Use `--allow-overwrite submission_metadata`, not `--force`. For a large database, run it in slices with `--start-date` and `--end-date`.

The backfill also records the inbox of each submission that has none (#696). It lists the inboxes of the submission's LE, and it records an inbox only if exactly one of them holds the submission.

### 5. Recompute the QC verdicts (#656)

Under GRZ_QC_Workflow 4.0.0, a metric fails only if the computed value is below the BfArM threshold. A deviation of more than 10% between the computed and the LE-provided value is reported, but does not fail the QC. The deviation counts in both directions, relative to the computed value.

New submissions get the new verdicts from the workflow. For the verdicts already stored, run `recompute-qc`. It reads the stored QC results and the `submission_metadata` from step 4, so it works for results of any workflow version:

```bash
grzctl db submission recompute-qc --year 2026 --quarter 3            # preview
grzctl db submission recompute-qc --year 2026 --quarter 3 --apply
```

Verdicts change in both directions. A failure caused only by a deviation becomes passed. A computed value below a threshold, which older workflows did not check, becomes failed. A rerun is safe, so run it again for the current quarter just before you generate the quarterly report.

### 6. Update your scripts

| Before                                                                | Now                                                                  |
| --------------------------------------------------------------------- | -------------------------------------------------------------------- |
| `grzctl <command> --config-file FILE`                                 | `grzctl --config FILE <command>`                                     |
| `grzctl download`, `grzctl clean`, `grzctl decrypt`                   | add `--inbox NAME`                                                   |
| `grzctl list`, `grzctl db sync-from-inbox`                            | add `--submitter-id LE_ID --inbox NAME`                              |
| `grzctl validate`                                                     | add `--submitter-id LE_ID`                                           |
| `grzctl submit`, `grzctl upload`                                      | removed, use grz-cli                                                 |
| `grzctl db submission populate --submission_date`                     | `--submission-date`                                                  |
| `db submission modify ID pseudonym VALUE`, `--ignore-field pseudonym` | `local_case_id` in place of `pseudonym`                              |
| `pseudonym` in JSON output                                            | now the case's psn; the old value is under `local_case_id`           |
| `grzctl dump-config` output                                           | YAML with masked secrets; `--reveal-secrets` shows them              |
| `db submission update --failure-reason network_error`, `upload_error` | `transfer_error`                                                     |

---

## ⚠️ Breaking / behavior changes

### grzctl: `db submission populate` (#681)

- Without `--submission-date`, the upload date is the S3 `LastModified` of `metadata/metadata.json` in the inbox. Pass `--inbox` if the submitter has several inboxes and the database has no inbox for the submission (#696).
- After `grzctl clean`, pass `--submission-date`.
- Replacing or removing a stored value needs `--force` or `--allow-overwrite FIELD`.

### grzctl: `db backfill` (#676, #677, #681)

- A `metadata.json` found in both archives counts as an error.
- A submission with an overwrite that no option allows is skipped whole.
- `--allow-overwrite donors` allows donor changes. Replacing a case link needs `--force`.

### grzctl: other commands

- `validate` fails an initial submission whose case already has a QC-passed initial submission, with the failure reason `duplicate_initial` (#633).
- `encrypt` and `archive` pick the archive from the consent in the metadata (#659).
- `db.known_public_keys` takes a list of keys, and a path to a file moves to `db.known_public_keys_file` (#693). The file may contain blank and `#` comment lines (#676).

### grzctl: quarterly report (#656)

- `3-Detailprüfung_*.tsv` has two more columns at the end: `detailed_qc_passed` per submission (`yes`, `no` or `not_performed`) and `deviation_exceeds_tolerance` per lab datum. With `--with-submission-ids`, `submission_id` moves two columns to the right.
- A submission that passes with a deviation stays in the table, with `detailed_qc_passed` = `yes`. `number_of_failed_qcs` counts only failures, not deviations.

### grzctl: failure reasons (#690)

- A failed step records a failure reason that names the cause and who has to act, see [Error handling](../../packages/grzctl/docs/error-handling.md). The new reasons are `detailed_qc_error`, `transfer_error`, `configuration_error`, `pruefbericht_generation_error`, `pruefbericht_rejected`, `submission_cleaned` and `interrupted`.
- `network_error` and `upload_error` are retired. Older states keep them, but `--failure-reason` no longer accepts them.
- SIGTERM stops grzctl the way Ctrl-C does, and the running step records `interrupted`.
- A rerun of `archive` for an archived submission records `ARCHIVED` instead of failing.

### grzctl: crypt4gh keys (#694)

- `decrypt` decrypts with the key of the submission's inbox, which the next section describes. So the inboxes of one LE may use different keys.
- Every key path must name a regular file. Otherwise every grzctl command stops.
- An inbox's `private_key_passphrase` comes before `C4GH_PASSPHRASE`.
- `archives.*.public_key` takes the public key inline, as `keys.grz_public_key` did in v4.0.0.

### grzctl: the inbox of a submission (#696)

- The database records the inbox that a submission came from. `download` records it, and so do `db sync-from-inbox` and `db submission populate --inbox`. `db backfill` records it for older submissions, see step 4.
- `download`, `clean` and `decrypt` take the inbox from `--inbox`, else from the database, else the LE's only inbox. `download` then also searches the LE's inboxes for the submission.
- These three commands read the database only with `--update-db`, the default. Then the database must be reachable.
- If no inbox resolves, `decrypt` fails and records `configuration_error`.
- `encrypt` signs the files that it re-encrypts for an archive with the private key of the submission's inbox. It takes the inbox from the database, else the LE's only inbox. If no inbox resolves, it signs with a random key and logs a warning.
- `decrypt --archive consented|non-consented` decrypts an archived submission. It takes the key from `--private-key-path`, else from the new optional `archives.<name>.private_key` or `private_key_path`. `--private-key-path` alone replaces the inbox key. Pass `--no-update-db`, so that the decrypt leaves the state of the submission unchanged.
- `list` and `db sync-from-inbox` need `--inbox` only if the LE has several inboxes.
- `db submission show` lists the inbox.

### grzctl: Prüfbericht config (#697)

- `pruefbericht.authorization_url`, `client_id`, `client_secret` and `api_base_url` are required.
- A config without one of them stops every grzctl command.

### grz-cli (#691, #694)

- grz-cli signs the encrypted files with the submitter private key if `keys.submitter_private_key` or `keys.submitter_private_key_path` is set. Otherwise it signs them with a random key, as before.
- LEs need not change anything.
- Decryption does not check the sender key, so grzctl decrypts these files as before.
- grz-cli takes the passphrase of the submitter private key from `keys.submitter_private_key_passphrase`, else from `C4GH_PASSPHRASE`, else from a prompt.

### Python API

- grz-common, grzctl: the secret config fields are `SecretStr | None`. Read them with `grz_common.models.base.get_secret_value()` (#680).
- grz-pydantic-models: `StrictIgnoringBaseModel` is renamed to `LosslessBaseModel` (#654).
- grz-db: the `withhold_destructive` and `has_pending_destructive` methods are removed. Use `SubmissionChangeSet.undeclared_destructive_changes()` (#681).
- grz-common: the expected exceptions derive from `grz_common.exceptions.GrzError`. `SubmissionValidationError` moves there from `grz_common.workers.submission`, and `grz_common.workers.download.DownloadError` is removed (#690).
- grz-common: `Crypt4GH.prepare_c4gh_keys`, `Crypt4GH.decrypt_file`, `Submission.encrypt`, `EncryptedSubmission.decrypt`, `Worker.encrypt` and `Worker.decrypt` take `X25519PrivateKey` and `X25519PublicKey` objects instead of paths. The `Crypt4GH` key loaders return these objects (#694).
- grz-common, grz-cli: `grz_common.models.keys` is removed. Import `KeyModel` and `KeyConfigModel` from `grz_cli.models.config`. The public key check is the type `grz_common.models.base.Crypt4GHPublicKey` (#694).
