# grz-gatekeeper

An API server for the GRZ submission process.

## TODO
- keep track of UploadId to be able to cancel and resume uploads (and get rid of orphaned files)

## Overview

Example GRZ side config:
```yaml
# gatekeeper.yaml
identifiers:
  grz: "GRZM00123"                                     # from the WES_tumor+germline example
  le: "123456789"                                      # from the WES_tumor+germline example

keys:
  grz_public_key_path: '/path/to/inbox.pub'
  grz_private_key_path: '/path/to/inbox.sec'
s3:
  api_base_url: "http://127.0.0.1:54321"                # default
  endpoint_url: "http://localhost:9000"                 # default in grz-watchdog testbed
  access_key: "Gd4G0VqhHTmi28twgsmF"                    # for the GRZ S3 user responsible for inbox access
  secret: "KEx8ABzui36NhyXbiQwYXxps3HahqD7EKL54v65f"    # for the GRZ S3 user responsible for inbox access
  bucket: test1                                         # default in grz-watchdog testbed

db:
  database_url: "sqlite:////path/to/submissions.sqlite" # e.g. from grz-watchdog testbed
  known_public_keys: "config/configs/known_keys"        # e.g. from grz-watchdog testbed
  author:
    name: "Alice"                                       # default in grz-watchdog testbed
    # THIS IS A DUMMY TEST KEY!
    private_key: |
      -----BEGIN OPENSSH PRIVATE KEY-----
      b3BlbnNzaC1rZXktdjEAAAAACmFlczI1Ni1jdHIAAAAGYmNyeXB0AAAAGAAAABBF9yxkGm
      ThDPJqANnAncZAAAAAGAAAAAEAAAAzAAAAC3NzaC1lZDI1NTE5AAAAIDu+8yDmRM755qm0
      A4oC9QALzyDgo4rIaOwF01p+ou8GAAAAkOWnOYjxMIsLYQdDzYJ0I2sbkXa2ppSzu+/dvh
      nkrMAJ9DAvvZk/7tykTcSxynBcf/DDLtb1mJTtuGaQ2MYtvxU0kOVB2sjWqiuNySO8FAKc
      NXGNF+XE/VVUeYmOGH5fGJaiIl32GUmhcOrz6Gb8WBSnB9CvH3FA5nmzSEqShU/vxd4UA7
      z5kwx1VgIvFS9XiQ==
      -----END OPENSSH PRIVATE KEY-----

    private_key_passphrase: "test"                      # default in grz-watchdog testbed

# new section for managing authentication
auth:
  secret_key: "f06adaba52642bd9178619fe8f615caf6a152a41b9061caa1ff022c8868185f1"  # used for jwt encoding
  algorithm: "HS256"
  access_token_expire_minutes: 120
  users:
    le-123456789:
      hashed_password: "$2b$12$sxNbW88d2o0aPK7kAMcGl.9dr1xlL57XNg2mWBV3YTOxFBNx9lzOm"  # as produced by `passlib.context.CryptContext(schemes=["bcrypt"], deprecated="auto").hash($SECRET)`
      disabled: false
```

Run on the GRZ side
```sh
grz-gatekeeper run --config-file gatekeeper.yaml
```

Example LE side config:
```yaml
# GRZM00123.config.yaml
identifiers:
  grz: "GRZM00123"                                     # from the WES_tumor+germline example
  le: "123456789"                                      # from the WES_tumor+germline example

keys:
  grz_public_key_path: '/path/to/inbox.pub'
s3:
  api_base_url: "http://127.0.0.1:54321"
  endpoint_url: "http://127.0.0.1:54321"  # unused if api_base_url is defined
  access_key: "le-123456789"              # re-use access_key as username for oauth2
  secret: "a_very_secret_password"        # re-use secret as password for oauth2; use envvars instead!
  bucket: test1                           # unused if api_base_url is defined
```

Upload to the GRZ via grz-cli
```sh
grz-cli upload --config-file GRZM00123.config.yaml --submission-dir /path/to/submission/dir
```