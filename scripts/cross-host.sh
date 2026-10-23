#!/usr/bin/env bash
#
# Check that the reference reproduces the tracked fixtures byte for byte on glibc: both codes take
# their maths from CORE-MATH, not the host libm, so the fixtures are host-independent (CLAUDE.md
# Rule 1, xtask/src/cross_host.rs). Exits non-zero, naming the files, if any record differs.
#
#   scripts/cross-host.sh
#
# Builds the Linux image in scripts/cross-host/, copies the tree into target/cross-host-tree/ (the
# container never sees the live checkout, so nothing in it can be overwritten; the copy keeps its own
# build cache between runs, and `cargo xtask` rebuilds the reference there whenever its inputs have
# changed) and runs `cargo xtask cross-host`.
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
tree="$repo_root/target/cross-host-tree"
image="yfoil-cross-host:ubuntu24"

command -v docker >/dev/null 2>&1 || { echo "error: docker is needed for the Linux side" >&2; exit 1; }
docker build -q -t "$image" "$repo_root/scripts/cross-host" >/dev/null
mkdir -p "$tree"
rsync -a --delete --exclude target --exclude .tmp --exclude .git --exclude 'studies/*/runs' "$repo_root/" "$tree/"
docker run --rm -v "$tree":/work -w /work -e CARGO_TARGET_DIR=/work/target "$image" \
    bash -c 'cargo run -q -p xtask -- cross-host'
