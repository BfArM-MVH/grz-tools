# End-to-end tests

These tests take a submission from the LE's upload to the GRZ's last step, with the real command line tools in their own processes.
The LE uploads with `grz-cli`, and the GRZ processes the submission in each of its two flows:

- `grzctl process`,
- the single commands from `download` to `clean`, with detailed QC by hand (`db should-qc`, `db submission update`), as GDCs run them without `process`.

The test then checks database, buckets, Prüfberichte and QC output.

```bash
uv run tox -e e2e                           # all scenarios
uv run tox -e e2e -- -k duplicate -n0       # some, in one process
```

They need docker or podman; `GRZ_E2E_CONTAINER_ENGINE` picks one.
Every scenario runs on SQLite and on PostgreSQL; the PostgreSQL half needs `pg_config` on the `PATH`, as the other database tests do.
The default `pytest` run skips them through the `e2e` marker.

## What runs where

| Part        | Runs as |
| ----------- | ------- |
| S3          | [SeaweedFS](https://github.com/seaweedfs/seaweedfs) in a container, one per test worker; each set of credentials has a user whose bucket policies grant exactly the [S3 permissions](../../packages/grzctl/README.md#s3-permissions) of the grzctl README. |
| Database    | `grz_db.testing`'s `migrated_db_connection`, a fresh one per test. |
| BfArM API   | An HTTPS mock in the test process, which a scenario can make reject the Prüfbericht. |
| Detailed QC | A shell command that lists the files it got, or fails. |
| Test data   | The submissions of [grz-mini-test-data](https://github.com/twrightsman/grz-mini-test-data), pinned by SHA-256 and cached in `.pytest_cache`. |
| Keys        | The mock keys in `tests/mock_files`. |

Mini files fit into one part, so the configs lower the part size to the 5 MiB that S3 allows at least, and two scenarios grow FASTQs beyond it.
SeaweedFS lets `s3:PutObject` abort a multipart upload, so the tests cannot show that `s3:AbortMultipartUpload` is needed; AWS and Ceph RGW keep the two apart.
They do check that no multipart upload is left open.

## Files

- `scenarios.py`: the scenarios.
  Each names a mini submission, how the LE changes and uploads it, the runs of the GRZ with the step where a run fails, and the outcome.
  `expected_states` derives the state history of each flow from that.
  A new scenario is a new entry in `SCENARIOS`, and runs in both flows.
- `test_scenarios.py`: runs a scenario and checks its outcome.
  For every scenario it also checks that the interrogation bucket keeps nothing, that the archive decrypts to the checksums in `metadata.json` (with the upstream crypt4gh library), that no decrypted copy is left unless the submission was selected for QC, and that grzctl logs no unexpected state transition.
- `harness.py`: the LE's and the GRZ's side, with their configs and commands.
- `services.py`: SeaweedFS with the permissions of each set of credentials, and the BfArM mock.
- `test_permissions_doc.py`: checks that the grzctl README lists the permissions that `services.py` grants; it runs without containers.
- `data.py`: download of the mini submissions.
