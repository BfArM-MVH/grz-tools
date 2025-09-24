# receive intermediate image name via build-arg
ARG SNAKEMAKE_BASE_IMAGE
# use intermediate image as foundation
# inherits all files, directories, and metadata (labels) (hopefully)
FROM ${SNAKEMAKE_BASE_IMAGE}

RUN apt-get update && \
    apt-get install -y --no-install-recommends curl && \
    rm -rf /var/lib/apt/lists/*

# install pixi into image
ENV PIXI_HOME="/root/.pixi"
ENV PIXI_VERSION="0.55.0"
RUN curl -fsSL https://pixi.sh/install.sh | bash

ENV PATH="${PIXI_HOME}/bin:${PATH}"

# install pixi-based runtime environment
# provides snakemake and global dependencies
WORKDIR /app
COPY pixi.toml pixi.lock ./
RUN pixi install

# copy workflow and configuration files
COPY workflow/ workflow/
COPY config/ workdir/config/

# set final entrypoint, use snakemake from pixi env
# overrides base image entrypoint
ENTRYPOINT ["pixi", "run", "snakemake", "--sdm", "conda", "--conda-prefix", "/conda-envs", "--directory", "/app/workdir"]

CMD []
