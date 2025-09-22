FROM ghcr.io/prefix-dev/pixi:bullseye

RUN apt-get update && \
    apt-get install -y ca-certificates curl && \
    rm -rf /var/lib/apt/lists/*
COPY local-certs/cert.pem /usr/local/share/ca-certificates/bfarm-mock-api.crt
RUN update-ca-certificates

WORKDIR /app
COPY pixi.toml pixi.toml
COPY pixi.lock pixi.lock
SHELL ["/bin/bash", "-c"]
RUN pixi install
COPY workflow/ workflow/
RUN mkdir /app/workdir
COPY config workdir/config
ENTRYPOINT ["pixi", "run", "snakemake", "--printshellcmds", "--cores", "2", "--sdm", "conda", "--directory", "/app/workdir"]