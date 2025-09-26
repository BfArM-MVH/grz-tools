#!/bin/bash
set -euo pipefail

# configuration
# final, usable image for grz-watchdog workflow
GRZ_WATCHDOG_WORKFLOW_IMAGE="grz-watchdog:latest"
# temporary image, contains conda envs and labels from snakemake
SNAKEMAKE_CONTAINERIZED_ENVS_IMAGE="grz-watchdog-conda-envs:latest"

# generate dockerfile from snakemake
snakemake --containerize > Dockerfile.snake.tmp

# check if an old Dockerfile exists and if it's identical to the new one
if [ -f Dockerfile.snake ] && cmp -s Dockerfile.snake Dockerfile.snake.tmp; then
  rm Dockerfile.snake.tmp
else
  mv Dockerfile.snake.tmp Dockerfile.snake
fi



# build intermediate image
podman build -t "${SNAKEMAKE_CONTAINERIZED_ENVS_IMAGE}" -f Dockerfile.snake .

# build final image, inherits from intermediate image
# pass intermediate image name as build-arg
podman build \
	--build-arg SNAKEMAKE_BASE_IMAGE="${SNAKEMAKE_CONTAINERIZED_ENVS_IMAGE}" \
	-t "${GRZ_WATCHDOG_WORKFLOW_IMAGE}" \
	-f Dockerfile .

# clean up intermediate artifacts
podman rmi "${SNAKEMAKE_CONTAINERIZED_ENVS_IMAGE}"
