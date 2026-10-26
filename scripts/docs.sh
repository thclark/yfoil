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
# Bump it deliberately, and check the rendered site before committing the bump. The pin moved
# from 0.0.60 to 0.0.71 (2026-10-23) for the exclude plugin (from 0.0.68), which keeps the
# generated table fragments (zensical.toml) out of the build; no page's article content changed.
set -euo pipefail

zensical_version="0.0.71"

# The ironlab-web bundle that draws the <ironlab-figure> embeds (docs/development/conventions/figures.md).
# Pinned by version and checksum, fetched from the ironlab release into a versioned directory
# (gitignored), so an upgrade changes every URL and no browser pairs a page with a stale engine.
ironlab_version="0.16.0"
ironlab_sha256="ad1160b9052f555eb060390b3fbd04116448c10888f53eac69ecd8dbd2e67380"

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

if ! command -v uv >/dev/null 2>&1; then
    echo "error: uv is not installed. Install it from https://docs.astral.sh/uv/ and re-run." >&2
    exit 1
fi

ironlab_dir="docs/javascripts/ironlab-$ironlab_version"
if [[ ! -f "$ironlab_dir/ironlab.js" ]]; then
    archive="ironlab-web-$ironlab_version.tar.gz"
    mkdir -p .tmp "$ironlab_dir"
    curl -fsSL -o ".tmp/$archive" \
        "https://github.com/thclark/ironlab/releases/download/$ironlab_version/$archive"
    echo "$ironlab_sha256  .tmp/$archive" | shasum -a 256 -c --quiet -
    tar xzf ".tmp/$archive" -C "$ironlab_dir"
fi

# The sidebar rules (icons, well-formed front matter) are checked before every build and serve:
# scripts/docs-check.py names each page that breaks one.
uv run --no-project --quiet python scripts/docs-check.py

uv run --no-project --quiet \
    --with "zensical==$zensical_version" \
    zensical "$@"

# Zensical does not fingerprint the site's own stylesheets, scripts and images, so a
# browser could pair a freshly deployed page with ones it cached from the last deploy.
# Stamp every reference to them with a content hash (see scripts/docs-cachebust.py).
if [[ "${1:-}" == "build" ]]; then
    uv run --no-project --quiet python scripts/docs-cachebust.py site
fi
