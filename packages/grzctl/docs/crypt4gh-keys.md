# Crypt4GH keys

Every encrypted file in a submission is a Crypt4GH file.
The GRZ config names the Crypt4GH keys that the different steps need:

- decrypting a downloaded submission needs the key of the inbox that it came from,
- re-encrypting a validated submission for the archive needs the archive's public key,
- signing the re-encrypted files needs the GRZ signing key.

This page explains where each key goes in the config, which key layout the primary
setup uses, and how to deviate from it.

## Key roles

| Key | Owned by | Where it appears | Used for |
| --- | -------- | ---------------- | -------- |
| GRZ key pair | GRZ | `leistungserbringer.<le>.inbox_buckets.<inbox>.private_key[_path]`, and `archives.signing_key[_path]` for signing | decrypting a submission, signing files for the archive |
| LE key pair | LE (submitter) | `keys.submitter_private_key_path` in the grz-cli config, `keys.grz_public_key` is the GRZ public key | signing for the inbox |
| Archive key pairs | GRZ | `archives.consented.public_key[_path]`, `archives.non_consented.public_key[_path]` | re-encrypting files for the archive |

A _key location_ is one of `private_key`, `private_key_path`, `public_key`, or `public_key_path`.
Every key location is validated the same way:

- `key` and `key_path` are mutual: setting both is a configuration error.
- A `key_path` must name an existing file.
- A required key location must set exactly one of the two fields.

A private key may be encrypted.
The passphrase comes from the first of:

1. the `..._passphrase` config field of the key,
2. the `C4GH_PASSPHRASE` environment variable,
3. an interactive prompt.

The LE's private key is a Crypt4GH key in the OpenSSH format.
`grz-cli encrypt` reads it from `keys.submitter_private_key_path` in the grz-cli config.
A passphrase-protected LE key uses the same `C4GH_PASSPHRASE` chain.

## Primary scenario: one GRZ key pair

The simplest layout uses a single GRZ key pair for every inbox and for signing.
The LE encrypts each submission to the GRZ public key.
The GRZ (grzctl) decrypts it with the matching private key, and re-encrypts it for the
archive with `archives.<consented|non_consented>.public_key`.
The GRZ private key signs all files that grzctl writes to an archive.
A YAML anchor defines the key once and reuses it everywhere.

```yaml
x-grz-key: &grz_key_path /etc/grzctl/keys/grz.sec

leistungserbringer:
  "260914050": # LE ID, the key of the leistungserbringer mapping
    alias: "GZD_LE" # optional
    inbox_buckets:
      inbox: # the name that --inbox and --inbox-name take
        endpoint_url: https://s3.example.org
        bucket: le-260914050 # optional, defaults to the inbox name
        private_key_path: *grz_key_path
  # …more LEs use the same anchor

archives:
  consented:
    s3:
      endpoint_url: https://s3.example.org
      bucket: grz-archive-consented
    public_key_path: /etc/grzctl/keys/archive-consented.pub
  non_consented:
    s3:
      endpoint_url: https://s3.example.org
      bucket: grz-archive-non-consented
    public_key_path: /etc/grzctl/keys/archive-non-consented.pub
  signing_key_path: *grz_key_path

db:
  database_url: sqlite:////var/lib/grzctl/grz.db
  author:
    name: alice
    private_key_path: /etc/grzctl/keys/author.sec

pruefbericht:
  base_url: https://example.org
  auth_url: https://example.org/auth
  client_id: …
  client_secret: …

identifiers:
  grz: GRZABC123
```

The archive public keys are the public halves of two separate key pairs.
The config holds only the public halves.
The private halves are not part of the grzctl config, so archived data is no longer open to the
LE that uploaded it.

## Variant: a separate signing key

A GRZ can sign the re-encrypted files with a key of its own that differs from the inbox key.
Set that key in `archives.signing_key` or `archives.signing_key_path` instead of reusing the
anchor:

```yaml
archives:
  …
  # signing_key: <inline crypt4gh private key>  # instead of the _path
  signing_key_path: /etc/grzctl/keys/signing.sec
```

The inbox private keys stay where the primary scenario puts them.
Only the signing role moves to the other key.

## Variant: several inboxes with different keys

Every inbox can name its own private key:

```yaml
leistungserbringer:
  "260914050":
    inbox_buckets:
      inbox: { endpoint_url: …, bucket: …, private_key_path: /etc/grzctl/keys/le-a.sec }
      inbox2: { endpoint_url: …, bucket: …, private_key_path: /etc/grzctl/keys/le-b.sec }
```

Decrypting picks the key by the inbox the submission came from.
When the database recorded that inbox, `decrypt` loads exactly its key.
Then each inbox of one LE may use its own key.
For a submission whose origin is not recorded, decryption falls back to the submitter.
Then the inboxes of one LE must all use one key.
Two inboxes use the same key when they name the same file or hold the same inline text,
for example through a YAML anchor.
