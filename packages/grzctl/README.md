# grzctl

Command-line tool for internal GRZ operations.

## Running a development version

1. Install [`uv`](https://docs.astral.sh/uv)
   - An easy way is to create a Conda environment containing `uv`.
2. Clone the `grz-tools` repository locally
3. From the repository root, use `uv run grzctl <grzctl options here>`
   - Alternatively, you can use `uv run --project path/to/repo grzctl <grzctl options here>` to run it from any directory.
     This is useful if your config uses relative paths and `grzctl` must therefore be run from a specific directory.

## Documentation

- [Registering change requests](docs/change-requests.md) — recording that a submission must be modified, deleted, or transferred, with the who/when/what audit trail.
- [Crypt4GH keys](docs/crypt4gh-keys.md): which Crypt4GH keys grz-cli and grzctl need, where they go in each config, and the key layouts of grzctl.
- [Example config](docs/example-config.yaml): a minimal grzctl config with placeholder values.
- [Processing a submission](docs/process.md): the `grzctl process` streaming pipeline, detailed QC prediction and selection, and recovery from partial failures.
- [Error handling](docs/error-handling.md): what a failed step records, the failure reasons and who acts on each.

## S3 permissions

Each set of credentials needs these actions on its buckets, for `grzctl process` and for the single commands from `download` to `clean` alike.
The names are those of bucket policies on AWS and Ceph RGW; other S3 servers may name them differently.
The end-to-end tests in `tests/e2e` run every flow with exactly these permissions.

| Credentials                                     | Bucket                | Actions |
| ----------------------------------------------- | --------------------- | ------- |
| `s3` of the LE's grz-cli config                 | inbox                 | `s3:GetObject`, `s3:PutObject`, `s3:AbortMultipartUpload`, `s3:ListBucket` |
| `leistungserbringer.<LE>.inbox_buckets.<inbox>` | inbox                 | `s3:GetObject`, `s3:PutObject`, `s3:DeleteObject`, `s3:ListBucket` |
| `archives.interrogation.s3`                     | interrogation         | `s3:GetObject`, `s3:PutObject`, `s3:DeleteObject`, `s3:AbortMultipartUpload`, `s3:ListBucket` |
| `archives.consented.s3`                         | consented archive     | `s3:GetObject`, `s3:PutObject`, `s3:AbortMultipartUpload`, `s3:ListBucket` |
| `archives.consented.s3`                         | interrogation         | `s3:GetObject`, `s3:ListBucket` |
| `archives.non_consented.s3`                     | non-consented archive | `s3:GetObject`, `s3:PutObject`, `s3:AbortMultipartUpload`, `s3:ListBucket` |
| `archives.non_consented.s3`                     | interrogation         | `s3:GetObject`, `s3:ListBucket` |

- `s3:ListBucket` lets S3 answer a missing object with 404 instead of 403.
  Without it, grz-cli and grzctl cannot tell a missing object from denied access, see [Error handling](docs/error-handling.md).
- The inbox credentials of the GRZ write and delete only in `grzctl clean` and in the clean step of `grzctl process`, which runs unless `--no-clean-inbox`.
  `grzctl inbox push-version` writes `version.json`.
  Credentials that run none of these need only `s3:GetObject` and `s3:ListBucket`.
- `grzctl process` copies each file from the interrogation bucket to its archive on the S3 server, with the credentials of that archive.
  So the interrogation bucket must be on the same endpoint as both archives, and both archives' credentials read it.
  The credentials of the interrogation bucket delete the staged copies afterwards.
- Write-only archive credentials, without `s3:GetObject` and `s3:ListBucket` on the archive, work too.
  `grzctl process` then asks the database whether an earlier run archived the submission, and `grzctl archive` warns and uploads.
- No command lists multipart uploads or their parts, so `s3:ListBucketMultipartUploads` and `s3:ListMultipartUploadParts` are not needed.

Grant `s3:AbortMultipartUpload`, not only the write actions; S3 servers often keep them apart.
A failed upload aborts its multipart upload, so that no parts stay behind.
Without the permission, the abort is skipped: the run still fails with the original error, but the uploaded parts stay in the bucket and cost storage.
A lifecycle rule that deletes incomplete multipart uploads after a few days covers that too.

## Sensitive configuration

Secret values (S3 secret keys, session tokens, key passphrases, Prüfbericht client secret) are read from the environment and stored as Pydantic `SecretStr`.
They are never written back to YAML, and `grzctl dump-config` prints `**********` in their place.
`grzctl dump-config --reveal-secrets` prints the plain values instead, so its YAML output loads back as a config file.

When a passphrase is not configured for a key, `grzctl process` falls back to the standard crypt4gh environment variable `C4GH_PASSPHRASE`.
If that is also unset, you will be prompted interactively.
Prefer an explicit `grz_...__private_key_passphrase` environment variable (or the `private_key_passphrase` config key) for automated runs.
