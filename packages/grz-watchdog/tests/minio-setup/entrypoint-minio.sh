#!/bin/sh
set -e

# marker file for healthcheck
SETUP_FLAG_FILE="/tmp/setup_complete"
rm -f "$SETUP_FLAG_FILE"

/usr/bin/minio server /data --console-address :9001 &
SERVER_PID=$!

# run setup script
/config/setup-minio.sh

# if successful, touch marker file
touch "$SETUP_FLAG_FILE"

wait $SERVER_PID
