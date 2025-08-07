# grz-gatekeeper

An API server for the GRZ submission process.

## TODO
- keep track of UploadId to be able to cancel and resume uploads (and get rid of orphaned files)

## Overview

Example config:
```yaml
identifiers:
  grz: "GRZM00123"                                     # from the WES_tumor+germline example
  le: "123456789"                                      # from the WES_tumor+germline example

keys:
  grz_public_key_path: '/path/to/inbox.pub'
  grz_private_key_path: '/path/to/inbox.sec'
s3:
  api_base_url: "http://127.0.0.1:54321"                # default
  endpoint_url: "http://localhost:9000"                 # default in grz-watchdog testbed
  access_key: "Gd4G0VqhHTmi28twgsmF"                    # default in grz-watchdog testbed
  secret: "KEx8ABzui36NhyXbiQwYXxps3HahqD7EKL54v65f"    # default in grz-watchdog testbed
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
```

Run on the GRZ side
```sh
grz-gatekeeper run --config-file gatekeeper.yaml
```

Upload to the GRZ via grz-cli
```sh
grz-cli upload --config-file inbox.yaml --submission-dir /path/to/submission/dir
```