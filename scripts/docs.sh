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
#
# The version is pinned here, and only here, so that a build is reproducible.
# Bump it deliberately, and check the rendered site before committing the bump.
set -euo pipefail

zensical_version="0.0.60"

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

if ! command -v uv >/dev/null 2>&1; then
    echo "error: uv is not installed. Install it from https://docs.astral.sh/uv/ and re-run." >&2
    exit 1
fi

uv run --no-project --quiet \
    --with "zensical==$zensical_version" \
    zensical "$@"

# Zensical does not fingerprint the site's own stylesheets, scripts and images, so a
# browser could pair a freshly deployed page with ones it cached from the last deploy.
# Stamp every reference to them with a content hash (see scripts/docs-cachebust.py).
if [[ "${1:-}" == "build" ]]; then
    uv run --no-project --quiet python scripts/docs-cachebust.py site
fi
