# The other host of the cross-host measurement (scripts/cross-host.sh): Ubuntu 24.04, glibc 2.39,
# gfortran and Rust, the same toolchain classes as the CI runner (.github/workflows/build-and-test.yml).
FROM ubuntu:24.04
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential gfortran libx11-dev libfontconfig1-dev pkg-config curl git ca-certificates python3 patch \
    && rm -rf /var/lib/apt/lists/*
ENV RUSTUP_HOME=/usr/local/rustup CARGO_HOME=/usr/local/cargo PATH=/usr/local/cargo/bin:$PATH
RUN curl -sSf https://sh.rustup.rs | sh -s -- -y --profile minimal --default-toolchain stable
WORKDIR /work
