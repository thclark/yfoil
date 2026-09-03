#!/usr/bin/env bash
# Build the XFOIL reference binaries from the vendored pristine source plus the
# tracked patch series. Never builds in-tree.
#
#   target/xfoil-ref/pristine/      pristine 6.99 + build-config patches only (double precision)
#   target/xfoil-ref/instrumented/  the above + instrumentation patches
#
# Usage: scripts/xfoil-build.sh [--verify] [--snan]
#   --verify   after building, run both binaries on the smoke case and assert
#              their numeric output is byte-identical (instrumentation is inert)
#   --snan     additionally build a pristine variant with -finit-real=snan
#              -ffpe-trap=invalid and run the smoke case (no uninitialised reads)
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SRC="$ROOT/third_party/xfoil-6.99"
PATCHES="$ROOT/xfoil-instrumentation"
OUT="$ROOT/target/xfoil-ref"
VERIFY=0; SNAN=0
for a in "$@"; do case "$a" in --verify) VERIFY=1;; --snan) SNAN=1;; *) echo "unknown arg $a" >&2; exit 2;; esac; done

# --- integrity of the vendored tree ------------------------------------------
( cd "$SRC" && shasum -a 256 -c --quiet "$ROOT/third_party/xfoil-6.99.sha256" ) \
  || { echo "third_party/xfoil-6.99 does not match its checksum manifest" >&2; exit 1; }

# --- stage a tree and apply a series -------------------------------------------
stage() { # stage <name> <extra-fflags> <series...>
  local name="$1"; shift; local extra="$1"; shift
  local dir="$OUT/$name"
  rm -rf "$dir"; mkdir -p "$OUT"; cp -R "$SRC" "$dir"
  for s in "$@"; do
    while IFS= read -r p; do
      [ -z "$p" ] && continue
      patch -s -p1 -d "$dir" < "$PATCHES/$p" || { echo "patch $p failed in $name" >&2; exit 1; }
    done < "$PATCHES/$s"
  done
  if [ -n "$extra" ]; then
    sed -i '' "s|^DBL = \(.*\)$|DBL = \1 $extra|" "$dir/bin/Makefile"
  fi
  # plotlib (double precision), then xfoil. The Makefile's final `cp ./xfoil xfoil`
  # install step fails with BINDIR=. even though the link succeeded, so success is
  # judged by the presence of a freshly linked binary, not by make's exit status.
  ( cd "$dir/plotlib" && { make -s clean >/dev/null 2>&1 || true; } && make -s > "$OUT/$name.plotlib.log" 2>&1 ) \
    || { echo "plotlib build failed for $name; see $OUT/$name.plotlib.log" >&2; exit 1; }
  ( cd "$dir/bin" && { make -s clean >/dev/null 2>&1 || true; } && rm -f xfoil && { make -s xfoil > "$OUT/$name.xfoil.log" 2>&1 || true; } && [ -x xfoil ] ) \
    || { echo "xfoil build failed for $name; see $OUT/$name.xfoil.log" >&2; exit 1; }
  echo "built $dir/bin/xfoil"
}

stage pristine     "" series.build
stage instrumented "" series.build series.instrument

# --- provenance manifest -------------------------------------------------------
{
  echo "source_sha256_manifest: $(shasum -a 256 "$ROOT/third_party/xfoil-6.99.sha256" | cut -d' ' -f1)"
  echo "series_build_sha256:    $(cat "$PATCHES"/build/*.patch | shasum -a 256 | cut -d' ' -f1)"
  echo "series_instr_sha256:    $(cat "$PATCHES"/instrument/*.patch | shasum -a 256 | cut -d' ' -f1)"
  echo "gfortran:               $(gfortran --version | head -1)"
  echo "fflags:                 $(grep -m1 '^FFLAGS' "$OUT/pristine/bin/Makefile") / $(grep -m1 '^DBL' "$OUT/pristine/bin/Makefile")"
  echo "host:                   $(uname -m)-$(uname -s)-$(uname -r)"
  echo "libm:                   $( (otool -L "$OUT/pristine/bin/xfoil" 2>/dev/null | grep -m1 libSystem | awk '{print $1, $NF}') || (ldd "$OUT/pristine/bin/xfoil" | grep -m1 libm) )"
} > "$OUT/manifest.txt"
cat "$OUT/manifest.txt"

# --- smoke case: identical geometry, alpha=2, Re=1e6, ITER 20 -----------------
smoke() { # smoke <bindir> <workdir>
  local bin="$1" wd="$2"; rm -rf "$wd"; mkdir -p "$wd"; cp "$ROOT/tests/fixtures/xfoil/smoke/panels.dat" "$wd/panels.dat"
  ( cd "$wd" && printf 'PLOP\nG F\n\nLOAD panels.dat\nOPER\nVISC 1000000\nITER 20\nALFA 2\nCPWR cp.dat\nDUMP bl.dat\n\nQUIT\n' | "$bin/xfoil" > stdout.txt 2>&1 || true )
}

if [ "$VERIFY" = 1 ]; then
  [ -f "$ROOT/tests/fixtures/xfoil/smoke/panels.dat" ] || { echo "smoke panels missing: tests/fixtures/xfoil/smoke/panels.dat" >&2; exit 1; }
  smoke "$OUT/pristine/bin" "$OUT/verify/pristine"
  smoke "$OUT/instrumented/bin" "$OUT/verify/instrumented"
  # Compare only the numeric products XFOIL writes in both builds (cp.dat, bl.dat, and the
  # OPER summary lines on stdout); instrumentation adds extra files and log lines by design.
  for f in cp.dat bl.dat; do
    cmp -s "$OUT/verify/pristine/$f" "$OUT/verify/instrumented/$f" \
      || { echo "INSTRUMENTATION NOT INERT: $f differs" >&2; exit 1; }
  done
  grep -E '^\s+a =|^\s+CL =|^\s+Cm =|^\s+CD =' "$OUT/verify/pristine/stdout.txt"     > "$OUT/verify/p.txt" || true
  grep -E '^\s+a =|^\s+CL =|^\s+Cm =|^\s+CD =' "$OUT/verify/instrumented/stdout.txt" > "$OUT/verify/i.txt" || true
  cmp -s "$OUT/verify/p.txt" "$OUT/verify/i.txt" || { echo "INSTRUMENTATION NOT INERT: OPER summary differs" >&2; exit 1; }
  echo "verify: instrumentation is inert on the smoke case (cp.dat, bl.dat, OPER summary byte-identical)"
fi

if [ "$SNAN" = 1 ]; then
  # -finit-real=snan poisons every uninitialised local REAL; -ffpe-trap=invalid aborts the
  # moment one is consumed. Success = the run completes with exit 0 and no backtrace.
  # gfortran's exit-time "Note: ... exceptions are signalling" lists raised-but-untrapped
  # flags (e.g. divide-by-zero) and is reported as information, not failure.
  stage snan "-g -fbacktrace -finit-real=snan -ffpe-trap=invalid" series.build
  rm -rf "$OUT/verify/snan"; mkdir -p "$OUT/verify/snan"; cp "$ROOT/tests/fixtures/xfoil/smoke/panels.dat" "$OUT/verify/snan/"
  set +e
  ( cd "$OUT/verify/snan" && printf 'PLOP\nG F\n\nLOAD panels.dat\nOPER\nVISC 1000000\nITER 20\nALFA 2\n\nQUIT\n' | "$OUT/snan/bin/xfoil" > stdout.txt 2>&1 ); rc=$?
  set -e
  if [ $rc -ne 0 ] || grep -q "Backtrace for this error" "$OUT/verify/snan/stdout.txt"; then
    echo "UNINITIALISED READ DETECTED (exit $rc) under -finit-real=snan -ffpe-trap=invalid:" >&2
    grep -A12 "Backtrace for this error" "$OUT/verify/snan/stdout.txt" >&2 || tail -20 "$OUT/verify/snan/stdout.txt" >&2
    exit 1
  fi
  echo "snan: no uninitialised REAL consumed on the smoke case (exit 0, no trap)"
  # Known and benign: IEEE_DIVIDE_BY_ZERO is raised in PLTINI (xplots.f:37, called from CPX
  # in OPER) — plot-scale setup that runs even with graphics off. Located by trapping
  # -ffpe-trap=zero and symbolising the backtrace; it is outside the solver.
  grep -i "exceptions are signalling" "$OUT/verify/snan/stdout.txt" | sed 's/^/snan info (PLTINI, benign): /' || true
fi
