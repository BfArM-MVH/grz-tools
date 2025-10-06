#!/bin/bash
set -euo pipefail

GRZ_WATCHDOG_WORKFLOW_IMAGE="grz-watchdog:latest"
# To be able to run snakemake --containerize, we have to ensure the env vars snakemake needs are at least present/empty
GRZ_KEYS__GRZ_PRIVATE_KEY_PASSPHRASE="" GRZ_S3__ACCESS_KEY="" GRZ_S3__SECRET="" \
  snakemake --containerize > Dockerfile.snake.tmp

if [ -f Dockerfile.snake ] && cmp -s Dockerfile.snake Dockerfile.snake.tmp; then
  rm Dockerfile.snake.tmp
else
  mv Dockerfile.snake.tmp Dockerfile.snake
fi

sed -i '1s/$/ AS conda-base/' Dockerfile.snake

sed -e '/# --- SNAKEMAKE-CONTAINERIZE-CONTENT-GOES-HERE ---/r Dockerfile.snake' \
    -e '/# --- SNAKEMAKE-CONTAINERIZE-CONTENT-GOES-HERE ---/d' \
    Dockerfile.template > Dockerfile

podman build -t "${GRZ_WATCHDOG_WORKFLOW_IMAGE}" -f Dockerfile .
