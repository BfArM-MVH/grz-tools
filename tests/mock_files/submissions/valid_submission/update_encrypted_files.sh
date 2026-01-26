#!/bin/bash
for i in files/*; do
    crypt4gh encrypt --recipient_pk "../../grz_mock_public_key.pub" < "$i" > "encrypted_files/$(basename $i).c4gh"
done
