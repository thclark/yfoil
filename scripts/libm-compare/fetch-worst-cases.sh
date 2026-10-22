#!/usr/bin/env bash
#
# Fetch CORE-MATH's worst-case (hard-to-round) argument lists for the twelve functions, at the
# commit third-party/core-math is vendored from, into target/libm-compare/worst-cases/ (118 MB,
# not tracked). The vendored subset has the sources only; the lists are what CORE-MATH's own
# `make check` runs.
#
#   scripts/libm-compare/fetch-worst-cases.sh [DIR]
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
out="${1:-$repo_root/target/libm-compare/worst-cases}"
commit="$(sed -n 's/.*commit=\([0-9a-f]*\).*/\1/p' "$repo_root/third-party/core-math.ORIGIN")"
[ -n "$commit" ] || { echo "error: no commit in third-party/core-math.ORIGIN" >&2; exit 1; }
base="https://gitlab.inria.fr/core-math/core-math/-/raw/$commit/src/binary64"

mkdir -p "$out"
for f in asin atan atan2 cos exp hypot log log10 pow sin sincos tanh; do
    curl -sSfL -o "$out/$f.wc" "$base/$f/$f.wc"
    printf '%-7s %8d lines\n' "$f" "$(wc -l <"$out/$f.wc")"
done
echo "$commit" >"$out/COMMIT"
echo "worst cases at core-math $commit in $out"
