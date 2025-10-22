# Genome Data Centers (GDC)

!!! warning "Only relevant for GDC staff"
    If you are staff at a clinic looking to upload submissions to a genome data center, see [these instructions](../clinic/index.md) instead.

## Configuration Files

Genome data centers are responsible for supplying the YAML configuration file for `grz-cli` to their associated clinics.

The purpose of this file is to provide values for the following:

- Crypt4GH public key for encryption
- Expected IDs for the clinic and data node (genome data center)
- S3 API parameters

A minimal configuration file might look like the following:

```yaml
keys:
  grz_public_key: |2
    -----BEGIN CRYPT4GH PUBLIC KEY-----
    i/TkbXjDDLnb0OvZ7VmF8GwRXXxTg1djpT8JC8GSfAw=
    -----END CRYPT4GH PUBLIC KEY-----

identifiers:
  grz: 'GRZABC123'
  le: '123456789'

s3:
  endpoint_url: 'https://your-s3-endpoint.com'
  bucket: 'grz-inbox-123'
```

Note that this does not include the required S3 `access_key` and `secret`.

There are two ways to provide these S3 secrets:

1. Integrated into the config file under the `s3` key (see [below](#s3)).
2. Provided separately and instructing clinics to define them in the usual AWS environment variables:
    - `AWS_ACCESS_KEY_ID`
    - `AWS_SECRET_ACCESS_KEY`

The first solution is the easiest but has the disadvantage of storing the secrets to disk, which may have security implications.

The following sections show documentation for the Pydantic models of each top-level section.
The attributes represent possible YAML keys.

### `keys`

::: grz_common.models.keys.KeyModel
    options:
      heading_level: 4
      show_root_heading: true
      show_bases: false

### `identifiers`

::: grz_common.models.identifiers.IdentifiersModel
    options:
      heading_level: 4
      show_root_heading: true
      show_bases: false

### `s3`

::: grz_common.models.s3.S3Options
    options:
      heading_level: 4
      show_root_heading: true
      show_bases: false
