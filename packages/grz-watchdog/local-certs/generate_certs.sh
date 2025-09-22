#!/bin/bash
set -euo pipefail

# ensure we can find openssl.cnf (right next to this script(
cd "$(dirname "$0")"

if [ -f "cert.pem" ]; then
    echo "Certificate 'cert.pem' already exists. Skipping."
    exit 0
fi

echo "Generating self-signed certificate with SAN…"

openssl req -x509 -newkey rsa:4096 -nodes \
  -out cert.pem -keyout key.pem -days 365 \
  -config openssl.cnf -extensions v3_req

echo "Certificate and key generated: cert.pem, key.pem"

