#!/usr/bin/env bash
#
# Reference NACA ordinates from the public-domain PDAS `naca456` program (Ralph L. Carmichael's
# revision of the NASA Langley programs of TM X-3284, TM X-3069 and TM-4741), for
# tests/fixtures/naca456/ and tests/naca456_series_tests.rs.
#
#   scripts/naca456-fixtures.sh            download (once), build the driver, write every case
#   scripts/naca456-fixtures.sh --tables   also regenerate src/geometry/series/six_series_tables.rs
#
# The distribution is fetched into target/naca456/ (never in-tree) and pinned by SHA-256; the
# driver (scripts/naca456-fixtures/driver.f90) links naca456's own modules unchanged and is built
# in double precision with fused multiply-add off, like the XFOIL reference. Cases come from
# scripts/naca456-fixtures/cases.txt; each becomes tests/fixtures/naca456/<slug>.dat at
# 17 significant figures, and manifest.json records the provenance.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="$ROOT/target/naca456"
SRC="$OUT/src"
ZIP_URL="https://www.pdas.com/packages/naca456.zip"
ZIP_SHA="327f9dae5552803990d0b683dd0e1d6c777562e4e3695b50d8cee80eb412f2c9"
FFLAGS="-fdefault-real-8 -ffp-contract=off -O0"
FIXTURES="$ROOT/tests/fixtures/naca456"
CASES="$ROOT/scripts/naca456-fixtures/cases.txt"
DRIVER="$ROOT/scripts/naca456-fixtures/driver.f90"

mkdir -p "$OUT" "$FIXTURES"
if [ ! -f "$OUT/naca456.zip" ]; then
  echo "downloading $ZIP_URL"
  curl -sSL -A "yfoil-fixtures" -o "$OUT/naca456.zip" "$ZIP_URL"
fi
got="$(shasum -a 256 "$OUT/naca456.zip" | cut -d' ' -f1)"
if [ "$got" != "$ZIP_SHA" ]; then
  echo "naca456.zip sha256 $got, expected $ZIP_SHA (the upstream package changed; review before re-pinning)" >&2
  exit 1
fi
rm -rf "$SRC"; mkdir -p "$SRC"
unzip -q -o "$OUT/naca456.zip" -d "$SRC" naca456.f90 nacax.f90 epspsi.f90 splprocs.f90 readme.txt input.txt

gfortran $FFLAGS -I "$SRC" -J "$OUT" -o "$OUT/driver" "$DRIVER"   # -J: module files stay out of the tree

if [ "${1:-}" = "--tables" ]; then
  python3 "$ROOT/scripts/naca456-fixtures/generate-tables.py" "$SRC/epspsi.f90" > "$ROOT/src/geometry/series/six_series_tables.rs"
  echo "wrote src/geometry/series/six_series_tables.rs"
fi

n=0
slugs=()
while IFS= read -r line; do
  case "$line" in ''|'#'*) continue;; esac
  slug="${line%% *}"
  body="${line#* }"
  printf '&NACA %s /\n' "$body" | "$OUT/driver" > "$FIXTURES/$slug.dat"
  slugs+=("\"$slug\"")
  n=$((n + 1))
done < "$CASES"

join() { local IFS=","; echo "$*"; }
cat > "$FIXTURES/manifest.json" <<JSON
{
  "reference": "PDAS naca456 (Carmichael 2001), NASA TM-4741 algorithm",
  "source_url": "$ZIP_URL",
  "source_sha256": "$ZIP_SHA",
  "driver": "scripts/naca456-fixtures/driver.f90",
  "driver_sha256": "$(shasum -a 256 "$DRIVER" | cut -d' ' -f1)",
  "cases_sha256": "$(shasum -a 256 "$CASES" | cut -d' ' -f1)",
  "gfortran": "$(gfortran --version | head -1)",
  "fflags": "$FFLAGS",
  "host": "$(uname -sm)",
  "generated_by": "scripts/naca456-fixtures.sh",
  "cases": [$(join "${slugs[@]}")]
}
JSON
echo "wrote $n fixtures to tests/fixtures/naca456/ (manifest.json)"
