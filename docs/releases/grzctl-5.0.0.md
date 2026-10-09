## Upgrade guide: grzctl 5.0.0

This guide is for GRZ operators. It lists what to do before and after the update, and what behaves differently. The changelog of each package lists every change, see the [releases](https://github.com/BfArM-MVH/grz-tools/releases).

The main changes:

- grzctl reads one config file for all commands (#635).
- grzctl tracks cases with `grzctl db case` (#633), see [Case tracking](../../packages/grzctl/docs/case-tracking.md).
- Detailed QC reports a deviation without failing (#656).
- A failed step records why it failed (#690), see [Error handling](../../packages/grzctl/docs/error-handling.md).
- The grzctl config names each crypt4gh key where grzctl uses it (#694), see [Crypt4GH keys](../../packages/grzctl/docs/crypt4gh-keys.md).
- The database records the inbox of each submission, so `download`, `clean` and `decrypt` find it without `--inbox` (#696).
- `grzctl inbox push-version` publishes the grz-cli version policy to every inbox (#667).

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

### 1. Update grzctl and GRZ_QC_Workflow

#### 1.1 grzctl v5.0.0

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

#### 1.2 GRZ_QC_Workflow v4.0.0

grzctl 4.0.0 rejects the reports of GRZ_QC_Workflow 4.0.0, so update grzctl first. Step 5 explains the new QC verdicts.

### 2. Write the unified config file (#635)

> ℹ️ This concerns the grzctl config only. The LE configs for grz-cli stay unchanged.

What changes:

- grzctl reads one file for all commands, by default `~/.config/grzctl/config.yaml`.
  grzctl no longer reads `~/.config/grz-cli/config.yaml`, the default of grzctl 4.0.0.
- `grzctl --config PATH <command>` replaces `--config-file` on each command.
  grzctl no longer merges several files.
- Every command checks the whole file.
  If one of the five required sections is missing, every command stops, also a command that does not use that section.
- Environment variables now override the file, see 2.3.

#### 2.1 Example

Merge your old config files into one file like this one.
The sections `leistungserbringer`, `archives`, `db`, `identifiers` and `pruefbericht` are required.
Every field without an `optional` comment is required.

```yaml
inbox_defaults: &inbox_defaults # optional; grzctl ignores this section, and the inboxes merge it with <<
  endpoint_url: https://s3.example.org
  private_key_path: /path/to/grz.sec # or private_key: the key inline
  # private_key_passphrase: ... # optional, else C4GH_PASSPHRASE, else a prompt
leistungserbringer: # at least one LE
  "123456789": # LE ID, quoted
    alias: "FOO" # optional
    inbox_buckets: # at least one inbox per LE
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
      access_key: ...
      secret: ...
    public_key_path: /path/to/consented.pub # or public_key: the key inline
    # private_key_path: /path/to/consented.sec # optional, only for decrypt --archive
  non_consented:
    s3:
      endpoint_url: https://s3.example.org
      bucket: grz-non-consented # must differ from the consented bucket
      access_key: ...
      secret: ...
    public_key_path: /path/to/non_consented.pub
db:
  database_url: postgresql+psycopg://...
  author:
    name: ... # no whitespace
    private_key_path: /path/to/author.sec # or private_key: the key inline
    # private_key_passphrase: ... # optional
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

The inbox names of one LE must differ in more than case, so not `Main` and `main`.

`access_key` and `secret` may be left out.
boto3 then looks for credentials itself, for example in `AWS_ACCESS_KEY_ID` and `AWS_SECRET_ACCESS_KEY`.
These credentials then apply to every inbox and archive without its own.

Every key field takes the key inline as `<name>` or a file as `<name>_path`.
Setting both is an error.
A private key has an optional `<name>_passphrase`.
[Crypt4GH keys](../../packages/grzctl/docs/crypt4gh-keys.md) lists every key and its field.

Check the result with `grzctl dump-config`.
It shows the values after the environment variables are applied.
Add `--reveal-secrets` to show the secrets (#680).

#### 2.2 Where the old settings go

| grzctl 4.0.0                                             | grzctl 5.0.0                                                                                             |
| -------------------------------------------------------- | -------------------------------------------------------------------------------------------------------- |
| one config file per command, merged with `--config-file` | one file, `--config PATH`                                                                                |
| `s3.*` of the `download`, `list` and `clean` configs     | `leistungserbringer.<LE_ID>.inbox_buckets.<inbox>.*`, one entry per inbox                                |
| `s3.bucket` of an inbox                                  | `leistungserbringer.<LE_ID>.inbox_buckets.<inbox>.bucket`, defaults to the inbox name                    |
| `s3.*` of the `archive` config                           | `archives.consented.s3.*` and `archives.non_consented.s3.*`                                              |
| `keys.grz_private_key_path`                              | `leistungserbringer.<LE_ID>.inbox_buckets.<inbox>.private_key_path` (#694, #696)                         |
| `C4GH_PASSPHRASE` for the GRZ private key                | `leistungserbringer.<LE_ID>.inbox_buckets.<inbox>.private_key_passphrase`; `C4GH_PASSPHRASE` still works |
| `keys.grz_public_key` / `keys.grz_public_key_path`       | `archives.<archive>.public_key` / `public_key_path`, one key per archive                                 |
| `keys.submitter_private_key_path`                        | removed; only grz-cli uses it                                                                            |
| `db.*`                                                   | unchanged, except `known_public_keys`                                                                    |
| `db.known_public_keys`, a path to a file                 | `db.known_public_keys_file`, or the keys as a list in `db.known_public_keys` (#693)                      |
| default `~/.config/grz-cli/known_public_keys`            | default `~/.config/grzctl/known_public_keys`; move the file                                              |
| `pruefbericht.*`, optional                               | `pruefbericht.*`, all four fields required (#697)                                                        |
| `identifiers.grz`, `identifiers.le`                      | `identifiers.grz`; `validate` takes the LE ID from `--submitter-id`                                      |

`grzctl encrypt` signs the files for the archives with the private key of the submission's inbox (#694, #696).

#### 2.3 Environment variables

Every field of the config file can also come from an environment variable.
In grzctl 4.0.0, the file won, and an environment variable only filled a field that the file lacked.
Now the environment variable wins.

The name is `GRZ_` and the field's path in the file, with `__` between the levels.
Upper or lower case does not matter.
A required field may be left out of the file if an environment variable sets it, for example to keep the secrets out of the file.

| Field in the file                                                          | Environment variable                                                             |
| -------------------------------------------------------------------------- | -------------------------------------------------------------------------------- |
| `leistungserbringer."123456789".inbox_buckets.main.access_key`             | `GRZ_LEISTUNGSERBRINGER__123456789__INBOX_BUCKETS__MAIN__ACCESS_KEY`             |
| `leistungserbringer."123456789".inbox_buckets.main.secret`                 | `GRZ_LEISTUNGSERBRINGER__123456789__INBOX_BUCKETS__MAIN__SECRET`                 |
| `leistungserbringer."123456789".inbox_buckets.main.private_key_path`       | `GRZ_LEISTUNGSERBRINGER__123456789__INBOX_BUCKETS__MAIN__PRIVATE_KEY_PATH`       |
| `leistungserbringer."123456789".inbox_buckets.main.private_key_passphrase` | `GRZ_LEISTUNGSERBRINGER__123456789__INBOX_BUCKETS__MAIN__PRIVATE_KEY_PASSPHRASE` |
| `archives.consented.s3.access_key`                                         | `GRZ_ARCHIVES__CONSENTED__S3__ACCESS_KEY`                                        |
| `archives.non_consented.s3.secret`                                         | `GRZ_ARCHIVES__NON_CONSENTED__S3__SECRET`                                        |
| `archives.consented.public_key`                                            | `GRZ_ARCHIVES__CONSENTED__PUBLIC_KEY`                                            |
| `db.database_url`                                                          | `GRZ_DB__DATABASE_URL`                                                           |
| `db.author.private_key_passphrase`                                         | `GRZ_DB__AUTHOR__PRIVATE_KEY_PASSPHRASE`                                         |
| `db.known_public_keys`                                                     | `GRZ_DB__KNOWN_PUBLIC_KEYS`, as a JSON list: `'["ssh-ed25519 AAAA... alice"]'`   |
| `pruefbericht.client_secret`                                               | `GRZ_PRUEFBERICHT__CLIENT_SECRET`                                                |
| `identifiers.grz`                                                          | `GRZ_IDENTIFIERS__GRZ`                                                           |

Rename the variables of grzctl 4.0.0.
grzctl ignores a name that matches no field, without a warning, so the old names do nothing now.

| grzctl 4.0.0                                                   | grzctl 5.0.0                                                                       |
| -------------------------------------------------------------- | ---------------------------------------------------------------------------------- |
| `GRZ_S3__*` of an inbox, for example `GRZ_S3__SECRET`          | `GRZ_LEISTUNGSERBRINGER__<LE_ID>__INBOX_BUCKETS__<INBOX>__*`                       |
| `GRZ_S3__*` of an archive                                      | `GRZ_ARCHIVES__CONSENTED__S3__*` / `GRZ_ARCHIVES__NON_CONSENTED__S3__*`            |
| `GRZ_KEYS__GRZ_PRIVATE_KEY_PATH`                               | `GRZ_LEISTUNGSERBRINGER__<LE_ID>__INBOX_BUCKETS__<INBOX>__PRIVATE_KEY_PATH`        |
| `GRZ_KEYS__GRZ_PUBLIC_KEY` / `GRZ_KEYS__GRZ_PUBLIC_KEY_PATH`   | `GRZ_ARCHIVES__<ARCHIVE>__PUBLIC_KEY` / `GRZ_ARCHIVES__<ARCHIVE>__PUBLIC_KEY_PATH` |
| `GRZ_DB__KNOWN_PUBLIC_KEYS`, a path                            | `GRZ_DB__KNOWN_PUBLIC_KEYS_FILE`                                                   |
| other `GRZ_DB__*`, `GRZ_PRUEFBERICHT__*`, `GRZ_IDENTIFIERS__*` | unchanged                                                                          |

Check the environment of your workers for `GRZ_*` variables.
A variable that the file shadowed in grzctl 4.0.0 now overrides the file.

Pitfalls:

- An environment variable cannot remove a field of the file.
  So to pass a key inline, for example with `..._PRIVATE_KEY`, remove `private_key_path` from the file.
  Otherwise grzctl stops with `Only one of private_key or private_key_path must be set`.
- Use the LE ID, not the alias.
- A variable whose path matches no LE or inbox of the file adds a new, incomplete entry, and every command stops.
- A shell accepts only letters, digits and `_` in a variable name, so `export` fails for an inbox named like `le-123456789`.
  Name the inbox `main` and set `bucket: le-123456789`.
  Or set the inbox as JSON one level up, which grzctl merges into the inbox of the file:
  `GRZ_LEISTUNGSERBRINGER__123456789__INBOX_BUCKETS='{"le-123456789": {"secret": "..."}}'`.
- An empty variable sets an empty value.
  Unset a variable instead.
- Environment variables only apply together with a config file.
  Without one, grzctl stops with `Missing required option '--config'`.

These environment variables are not config fields, and they stay as they are:
`C4GH_PASSPHRASE` for every private key without its own passphrase, `GRZ_PRUEFBERICHT_ACCESS_TOKEN` for `pruefbericht submit --token`, `GRZCTL_SHOULD_QC_SALT` and `GRZCTL_QC_WORKFLOW_VERSION`.

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

### 7. Publish the grz-cli version policy (#667)

grz-cli reads `version.json` from the inbox before `upload` and `submit`. It stops if its version is too old, or if the inbox has no `version.json`. grzctl 5.0.0 ships a policy file, and this command publishes it to every configured inbox:

```sh
grzctl inbox push-version
```

The command replaces the `version.json` in each inbox. The policy recommends grz-cli 3.0.0, and from 2026-10-19 on it requires grz-cli 2.0.0.

---

## ⚠️ Breaking / behavior changes

### grzctl: `db submission populate` (#681)

- Without `--submission-date`, `populate` keeps the stored upload date (#690). If the database holds none, the upload date is the S3 `LastModified` of `metadata/metadata.json` in the inbox. Pass `--inbox` if the submitter has several inboxes and the database has no inbox for the submission (#696).
- After `grzctl clean`, pass `--submission-date` for a submission without a stored upload date.
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
- An inbox's `private_key_passphrase` comes before `C4GH_PASSPHRASE`. Prefer the grzctl setting to `C4GH_PASSPHRASE`, for example the environment variable `GRZ_LEISTUNGSERBRINGER__123456789__INBOX_BUCKETS__MAIN__PRIVATE_KEY_PASSPHRASE`. `C4GH_PASSPHRASE` applies to every key that has no passphrase of its own.
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
