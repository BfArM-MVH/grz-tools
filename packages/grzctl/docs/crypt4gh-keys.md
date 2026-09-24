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
| GRZ private key | grzctl: `leistungserbringer.<LE ID>.inbox_buckets.<inbox>.private_key[_path]` | grzctl decrypts a submission from that inbox. |
| GRZ signing key | grzctl: `archives.signing_key[_path]` | grzctl signs the files that it re-encrypts for an archive. |
| Archive public keys | grzctl: `archives.consented.public_key[_path]`, `archives.non_consented.public_key[_path]` | grzctl re-encrypts a submission for the matching archive. |

`<name>[_path]` stands for the two fields `<name>` and `<name>_path`.
The GRZ gives its public key to its LEs.
In the primary setup, the GRZ signing key is the GRZ private key.
No grzctl command needs the private keys of the archives, so they are in no config.

`db.author.private_key[_path]` in the grzctl config is no Crypt4GH key.
It signs the submission states in the database, and this page does not cover it.

## Key fields

The field `<name>` holds the key inline, and the field `<name>_path` names a file with the key.
Setting both is a configuration error.
Every key except the LE private key is required.
For a required key, one of the two fields must be set.
A `_path` field must name an existing regular file.

An inline public key must be in the Crypt4GH format.
A public key file and a private key can also be in the OpenSSH format, as an ed25519 key.

A private key can have a passphrase.
grzctl takes the passphrase from the first of:

1. the key's `<name>_passphrase` field, for example `archives.signing_key_passphrase`,
2. the `C4GH_PASSPHRASE` environment variable,
3. a prompt.

grz-cli has no passphrase field, so it starts with `C4GH_PASSPHRASE`.

## grz-cli (LE)

```yaml
keys:
  grz_public_key_path: /path/to/grz.pub # or grz_public_key with the key inline
  submitter_private_key_path: /path/to/le.sec # optional, or submitter_private_key with the key inline
```

`grz-cli encrypt` encrypts every file for the GRZ public key.
It signs the files with the LE private key if one is set, and with a random key otherwise.

## grzctl (GRZ): one GRZ key pair

The primary setup uses one GRZ key pair for all inboxes of all LEs, and for signing.
A YAML anchor names the key file once, and the other fields reuse it.
The two archives have one key pair each, and the config holds only their public keys.
The example leaves out the S3 credentials.

```yaml
leistungserbringer:
  "123456789": # LE ID
    alias: "FOO" # optional
    inbox_buckets:
      inbox: # the name that --inbox takes
        endpoint_url: https://s3.example.org
        bucket: le-123456789 # optional, defaults to the inbox name
        private_key_path: &grz_key /path/to/grz.sec
  "000000000":
    inbox_buckets:
      inbox:
        endpoint_url: https://s3.example.org
        bucket: le-000000000
        private_key_path: *grz_key

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
  signing_key_path: *grz_key

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

If the GRZ private key has a passphrase, anchor `private_key_passphrase` and `signing_key_passphrase` the same way, or set `C4GH_PASSPHRASE`.

## grzctl (GRZ): a separate signing key

grzctl can sign with a key other than the inbox key.
Name that key in `archives.signing_key_path`, or put it inline in `archives.signing_key`, instead of the anchor:

```yaml
archives:
  ...
  signing_key_path: /path/to/signing.sec
```

## grzctl (GRZ): different keys per inbox

Each inbox can name its own private key.
grzctl decrypt picks the key by the submitter, because a submission's metadata names its submitter (LE) but not its inbox.
So all inboxes of one LE must use one key.
Different LEs can use different keys.
Two inboxes use one key when they name the same file or hold the same inline text, for example through a YAML anchor.
grzctl decrypt refuses to decrypt for an LE whose inboxes use more than one key.

```yaml
leistungserbringer:
  "123456789":
    inbox_buckets:
      inbox: { endpoint_url: ..., bucket: le-123456789, private_key_path: &key_a /path/to/key-a.sec }
      inbox2: { endpoint_url: ..., bucket: le-123456789-2, private_key_path: *key_a }
  "000000000":
    inbox_buckets:
      inbox: { endpoint_url: ..., bucket: le-000000000, private_key_path: /path/to/key-b.sec }
```
