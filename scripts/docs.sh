#!/usr/bin/env bash
#
# Build or preview the yFoil documentation site.
#
#   scripts/docs.sh serve          # live preview on http://localhost:8000
#   scripts/docs.sh serve --open   # ... and open a browser
#   scripts/docs.sh build --clean  # static site into site/
#
# Zensical is a Python package with a compiled Rust core. It is not installed
# into the repo: uv fetches the pinned version into its own cache on first run,
# so there is no virtualenv to manage and no Python state in the working tree.
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

if ! command -v uv >/dev/null 2>&1; then
    echo "error: uv is not installed. See https://docs.astral.sh/uv/" >&2
    echo "       (or install Zensical yourself: pip install -r requirements-docs.txt)" >&2
    exit 1
fi

exec uv run --no-project --quiet \
    --with-requirements requirements-docs.txt \
    zensical "$@"
