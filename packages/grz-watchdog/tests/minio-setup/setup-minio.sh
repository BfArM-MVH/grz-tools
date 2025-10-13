#!/bin/sh
set -e

until /usr/bin/mc alias set adm http://localhost:9000 "$MINIO_ROOT_USER" "$MINIO_ROOT_PASSWORD"; do
    echo "Waiting for Minio server..."
    sleep 1
done

/usr/bin/mc admin policy add adm readwrite /config/policy.readwrite.json
/usr/bin/mc admin user add adm "$MINIO_TEST_USER_ACCESS_KEY" "$MINIO_TEST_USER_SECRET_KEY"
/usr/bin/mc admin policy attach adm readwrite --user "$MINIO_TEST_USER_ACCESS_KEY"

exit 0
