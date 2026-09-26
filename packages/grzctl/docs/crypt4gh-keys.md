# Crypt4GH keys

The files of a submission are Crypt4GH files.
The LE encrypts them for the GRZ with grz-cli.
The GRZ decrypts them with grzctl and re-encrypts them for an archive.
This page lists the Crypt4GH keys that each tool needs, and where they go in its config.

## Which key goes where

| Key | Tool and config field | Used for |
| --- | --------------------- | -------- |
| GRZ public key | grz-cli: `keys.grz_public_key[_path]` | The LE encrypts a submission for the GRZ. |
| LE private key | grz-cli: `keys.submitter_private_key[_path]`, optional | grz-cli signs the files that it encrypts. |
| GRZ private key | grzctl: `leistungserbringer.<LE ID>.inbox_buckets.<inbox>.private_key[_path]` | grzctl decrypts a submission from that inbox. `grzctl encrypt` signs the files that it re-encrypts for an archive with this key, or with a random key if no inbox resolves. |
| Archive public keys | grzctl: `archives.consented.public_key[_path]`, `archives.non_consented.public_key[_path]` | grzctl re-encrypts a submission for the matching archive. |
| Archive private keys | grzctl: `archives.consented.private_key[_path]`, `archives.non_consented.private_key[_path]`, optional | `grzctl decrypt --archive` decrypts an archived submission. |

`<name>[_path]` stands for the two fields `<name>` and `<name>_path`.
The GRZ gives its public key to its LEs.
The archive private keys are optional. Only `grzctl decrypt --archive` uses them.
They decrypt the whole archive, so set them only on a host that needs to decrypt archived submissions.

`db.author.private_key[_path]` in the grzctl config is no Crypt4GH key.
It signs the submission states in the database, and this page does not cover it.

## Key fields

The field `<name>` holds the key inline, and the field `<name>_path` names a file with the key.
Setting both is a configuration error.
Every key except the LE private key and the archive private keys is required.
For a required key, one of the two fields must be set.
A `_path` field must name an existing regular file.

An inline public key must be in the Crypt4GH format.
A public key file and a private key can also be in the OpenSSH format, as an ed25519 key.

A private key can have a passphrase.
grzctl and grz-cli take the passphrase from the first of:

1. the key's `<name>_passphrase` field, for example `private_key_passphrase` of an inbox, or `keys.submitter_private_key_passphrase` in grz-cli,
2. the `C4GH_PASSPHRASE` environment variable,
3. a prompt.

## grz-cli (LE)

```yaml
keys:
  grz_public_key_path: /path/to/grz.pub # or grz_public_key with the key inline
  submitter_private_key_path: /path/to/le.sec # optional, or submitter_private_key with the key inline
  submitter_private_key_passphrase: ... # optional
```

`grz-cli encrypt` encrypts every file for the GRZ public key.
It signs the files with the LE private key if one is set, and with a random key otherwise.

## grzctl (GRZ): one GRZ key pair

The primary setup uses one GRZ key pair for all inboxes of all LEs.
The top-level mapping `inbox_defaults` holds the fields that all inboxes share, and each inbox merges it with `<<`.
The two archives have one key pair each, and the example config holds only their public keys.
The example leaves out the S3 credentials.

```yaml
inbox_defaults: &inbox_defaults # grzctl ignores this section; the inboxes merge it with <<
  endpoint_url: https://s3.example.org
  private_key_path: /path/to/grz.sec

leistungserbringer:
  "123456789": # LE ID
    alias: "FOO" # optional
    inbox_buckets:
      main: # the inbox name, which --inbox takes
        <<: *inbox_defaults
        bucket: le-123456789 # optional, defaults to the inbox name
  "000000000":
    inbox_buckets:
      main:
        <<: *inbox_defaults
        bucket: le-000000000

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

pruefbericht:
  authorization_url: https://...
  client_id: ...
  client_secret: ...
  api_base_url: https://...

identifiers:
  grz: GRZX00000
```

If the GRZ private key has a passphrase, put `private_key_passphrase` into `inbox_defaults` too, or set `C4GH_PASSPHRASE`.

## grzctl (GRZ): different keys per inbox

Each inbox can name its own private key.
This includes the inboxes of one LE.
Without `--archive` and `--private-key-path`, `grzctl decrypt` decrypts a submission with the key of the inbox that the submission came from.
grzctl looks up the inbox under the LE that the submission's metadata names.
It takes the inbox from `--inbox`.
Without `--inbox`, it takes the inbox recorded in the database, else the LE's only inbox.
`grzctl download` records the inbox, and `grzctl db backfill` records it for older submissions.

`grzctl encrypt` signs the files that it re-encrypts for an archive with the key of the same inbox.
It has no `--inbox`, so it takes the recorded inbox, else the LE's only inbox.
If neither resolves, it signs the files with a random key and logs a warning.

```yaml
leistungserbringer:
  "123456789":
    inbox_buckets:
      main: { endpoint_url: ..., bucket: le-123456789, private_key_path: /path/to/key-a.sec }
      inbox2: { endpoint_url: ..., bucket: le-123456789-2, private_key_path: /path/to/key-b.sec }
  "000000000":
    inbox_buckets:
      main: { endpoint_url: ..., bucket: le-000000000, private_key_path: /path/to/key-c.sec }
```

## grzctl (GRZ): decrypt an archived submission

An archived submission is encrypted for its archive.
`grzctl decrypt --archive` decrypts it with the archive's private key from the config:

```yaml
archives:
  consented:
    ...
    private_key_path: /path/to/consented.sec
```

```bash
grzctl decrypt --archive consented --no-update-db --submission-dir /path/to/submission
```

`--private-key-path` takes the key file on the command line instead, so the config needs no archive private key.
Its passphrase comes from `C4GH_PASSPHRASE`, else a prompt.

```bash
grzctl decrypt --private-key-path /path/to/consented.sec --no-update-db --submission-dir /path/to/submission
```

`--no-update-db` keeps the state of the submission in the database unchanged.
