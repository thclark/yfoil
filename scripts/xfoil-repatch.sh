#!/usr/bin/env bash
# Regenerate the instrumentation patch series from the instrumented working copy.
#
# Workflow for changing XFOIL instrumentation (CLAUDE.md Rule 3):
#   1. scripts/xfoil-build.sh                      # stages target/xfoil-ref/instrumented from pristine + patches
#   2. edit target/xfoil-ref/instrumented/src/*.f  # add/adjust WRITE statements (ES24.16, cwd-relative files)
#   3. scripts/xfoil-repatch.sh                    # rewrite xfoil-instrumentation/instrument/*.patch from the diff
#   4. scripts/xfoil-build.sh --verify             # rebuild both references, prove instrumentation still inert
#   5. cargo xtask fixtures                        # regenerate fixtures
#
# The build-config patches (xfoil-instrumentation/build/) are not touched here.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
P="$ROOT/third_party/xfoil-6.99"; W="$ROOT/target/xfoil-ref/instrumented"; OUT="$ROOT/xfoil-instrumentation/instrument"
[ -d "$W/src" ] || { echo "no instrumented working copy at $W; run scripts/xfoil-build.sh first" >&2; exit 1; }
# a build-config-only Makefile to diff the xlog hooks against
B="$(mktemp -d)"; cp -R "$P" "$B/tree"
while IFS= read -r p; do [ -n "$p" ] && patch -s -p1 -d "$B/tree" < "$ROOT/xfoil-instrumentation/$p"; done < "$ROOT/xfoil-instrumentation/series.build"
hdr() { sed -e "1s|.*|--- $1|" -e "2s|.*|+++ $2|"; }
diff -u "$B/tree/bin/Makefile" "$W/bin/Makefile" | hdr a/bin/Makefile b/bin/Makefile > "$OUT/10-makefile-xlog.patch" || true
{ diff -u /dev/null "$W/src/xlog.f" | hdr /dev/null b/src/xlog.f
  diff -u /dev/null "$W/src/XLOG.INC" | hdr /dev/null b/src/XLOG.INC
  diff -u "$P/src/XBL.INC" "$W/src/XBL.INC" | hdr a/src/XBL.INC b/src/XBL.INC; } > "$OUT/11-xlog-harness.patch" || true
n=12
for f in xfoil.f xgdes.f xpanel.f xbl.f xblsys.f xsolve.f xoper.f; do
  diff -u "$P/src/$f" "$W/src/$f" | hdr "a/src/$f" "b/src/$f" > "$OUT/$n-${f%.f}-dump.patch" || true
  n=$((n+1))
done
rm -rf "$B"
# guard: every numeric format on added lines must be ES24.16 (17 significant figures)
if grep -hoE "^\+.*[0-9]*(E|F|G)[0-9]+\.[0-9]+" "$OUT"/*.patch | grep -oE "[0-9]*(E|F|G)[0-9]+\.[0-9]+" | grep -vq "ES24.16"; then
  echo "WARNING: non-ES24.16 numeric formats on added lines:" >&2
  grep -hoE "^\+.*[0-9]*(E|F|G)[0-9]+\.[0-9]+" "$OUT"/*.patch | grep -oE "[0-9]*(E|F|G)[0-9]+\.[0-9]+" | grep -v "ES24.16" | sort | uniq -c >&2
fi
ls "$OUT"/*.patch | sort > "$ROOT/xfoil-instrumentation/series.instrument"
sed -i '' "s|^$ROOT/xfoil-instrumentation/||" "$ROOT/xfoil-instrumentation/series.instrument"
echo "regenerated $(wc -l < "$ROOT/xfoil-instrumentation/series.instrument") instrumentation patches"
