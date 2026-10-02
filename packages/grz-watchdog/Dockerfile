# Both stages use the same pixi, because the final stage runs the environment that the builder installs.
ARG PIXI_VERSION=0.75.0

FROM ghcr.io/prefix-dev/pixi:${PIXI_VERSION} AS builder

RUN apt-get update && apt-get install -y build-essential curl && rm -rf /var/lib/apt/lists/*

ENV RUSTUP_HOME=/usr/local/rustup CARGO_HOME=/usr/local/cargo PATH=/usr/local/cargo/bin:$PATH
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y

WORKDIR /build
COPY packages/grz-check/ ./
RUN cargo build --release --locked

WORKDIR /workspace/packages/grz-watchdog

COPY packages/grz-watchdog/pixi.toml packages/grz-watchdog/pixi.lock ./
COPY packages/grz-pydantic-models/ /workspace/packages/grz-pydantic-models/
COPY packages/grz-common/ /workspace/packages/grz-common/
COPY packages/grz-db/ /workspace/packages/grz-db/
COPY packages/grz-cli/ /workspace/packages/grz-cli/
COPY packages/grzctl/ /workspace/packages/grzctl/
COPY packages/grz-check/ /workspace/packages/grz-check/

# The setup task installs the grz packages into the environment. Building grz-check needs the Rust toolchain of this stage.
RUN pixi install --frozen && CARGO_TARGET_DIR=/build/target pixi run --frozen setup

FROM ghcr.io/prefix-dev/pixi:${PIXI_VERSION}

COPY --from=builder /build/target/release/grz-check /usr/local/bin/
COPY --from=builder /workspace /workspace

WORKDIR /workspace/packages/grz-watchdog

COPY packages/grz-watchdog/workflow /grz-watchdog/workflow
COPY packages/grz-watchdog/container_utils/scripts/ /usr/local/bin/
RUN chmod +x /usr/local/bin/*

RUN mkdir -p /workdir
WORKDIR /workdir

ENTRYPOINT ["/usr/local/bin/daemon"]
