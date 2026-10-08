#!/usr/bin/env bash
#
# Measure the reference's own spread between the fixtures' host and glibc, per tracked run case:
# the floor of the run tests' cross-host gate (CLAUDE.md Rule 1, xtask/src/cross_host.rs).
#
#   scripts/cross-host.sh        writes tests/fixtures/cross-host/spread.json
#
# Builds the Linux image in scripts/cross-host/, copies the tree into target/cross-host-tree/ (the
# container never sees the live checkout, so nothing in it can be overwritten; the copy keeps its own
# build cache between runs), builds the instrumented reference there and runs
# `cargo xtask cross-host`, then copies the measurement back. Rerun it whenever a tracked run case's
# fixture is regenerated: the tests refuse a measurement whose records no longer match.
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
tree="$repo_root/target/cross-host-tree"
image="yfoil-cross-host:ubuntu24"

command -v docker >/dev/null 2>&1 || { echo "error: docker is needed for the Linux side" >&2; exit 1; }
docker build -q -t "$image" "$repo_root/scripts/cross-host" >/dev/null
mkdir -p "$tree"
rsync -a --delete --exclude target --exclude .tmp --exclude .git --exclude 'scripts/*/runs' "$repo_root/" "$tree/"
docker run --rm -v "$tree":/work -w /work -e CARGO_TARGET_DIR=/work/target "$image" \
    bash -c 'cargo run -q -p xtask -- cross-host'
mkdir -p "$repo_root/tests/fixtures/cross-host"
cp "$tree/tests/fixtures/cross-host/spread.json" "$repo_root/tests/fixtures/cross-host/spread.json"
echo "tests/fixtures/cross-host/spread.json written"
