# grzctl

Command-line tool for internal GRZ operations.

## Running a development version

1. Install [`uv`](https://docs.astral.sh/uv)
  - An easy way is to create a Conda environment containing `uv`.
2. Clone the `grz-tools` repository locally
3. From the repository root, use `uv run grzctl <grzctl options here>`
  - Alternatively, you can use `uv run --project path/to/repo grzctl <grzctl options here>` to run it from any directory.
    This is useful if your config uses relative paths and `grzctl` must therefore be run from a specific directory.

## S3 permissions for `grzctl process`

`grzctl process` uploads to the archive buckets using multipart uploads. The
credentials you give it must be allowed to **abort** a multipart upload, not just
to write — these are often separate permissions.

This matters when an upload fails partway through: the tool tries to abort the
upload so no partial object is left behind. If the credentials lack abort
permission, the abort is skipped — processing still fails with the original
error, but the already-uploaded parts stay in the bucket and you can't remove
them yourself. Over time these orphaned parts pile up and cost storage.

Confusingly, S3 providers split multipart permissions in non-obvious ways (AWS,
for example, has both object-level and bucket-level multipart actions), and our
archive backend is Ceph/RGW rather than AWS, so check your provider's own
documentation for the exact action names. The practical point holds everywhere:
it's easy to grant write while leaving abort ungranted (see
[velero-io/velero#416](https://github.com/velero-io/velero/issues/416) for one
instance), so verify abort works — or rely on a bucket lifecycle rule that
deletes incomplete multipart uploads after a few days.

