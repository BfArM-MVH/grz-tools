FROM ghcr.io/prefix-dev/pixi:0.55.0-bullseye

RUN apt-get update && \
    apt-get install -y --no-install-recommends ca-certificates curl && \
    rm -rf /var/lib/apt/lists/*
COPY local-certs/cert.pem /usr/local/share/ca-certificates/bfarm-mock-api.crt
RUN update-ca-certificates

WORKDIR /app
COPY pixi.toml pixi.toml
COPY pixi.lock pixi.lock

RUN pixi install

COPY workflow/ workflow/
COPY config workdir/config

ENTRYPOINT ["pixi", "run", "snakemake", "--directory", "/app/workdir"]

CMD []
