#!/usr/bin/env bash
#
# Draw the branch-coverage study's occupancy figures (branch map, case matrix) from its run folder
# with matplotlib: the only figures ironlab cannot yet express (docs/contributing/conventions/figures.md; every
# other figure is an ironlab figure built in Rust).
#
#   scripts/figures/render.sh <study> <run_dir>      e.g. branch-coverage runs/2026-10-15T14-17-40Z
#
# The study's Rust program writes every number into the run folder; scripts/<study>/plot.py only
# presents them (scripts/figures/style.py holds the shared conventions). matplotlib is pinned here
# and fetched by uv into its own cache on first use, as scripts/docs.sh does for zensical; there is
# no Python environment in the working tree.
set -euo pipefail

matplotlib_version="3.11.1"

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
study="${1:?usage: render.sh <study> <run_dir>}"
run_dir="${2:?usage: render.sh <study> <run_dir>}"

if ! command -v uv >/dev/null 2>&1; then
    echo "error: uv is not installed. Install it from https://docs.astral.sh/uv/ and re-run." >&2
    exit 1
fi

cd "$repo_root"
exec uv run --no-project --quiet --with "matplotlib==$matplotlib_version" \
    python "scripts/$study/plot.py" "$run_dir"
