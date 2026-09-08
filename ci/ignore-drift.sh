#!/usr/bin/env bash
# Regenerate the list of #[ignore]d tests and fail if it differs from tests/IGNORED.txt.
# Every ignore must carry a reason naming its unblocking stage: #[ignore = "S7: ..."]
set -euo pipefail
cd "$(dirname "$0")/.."
tmp=$(mktemp)
grep -rnE '#\[ignore' tests src --include='*.rs' | sed 's/^\([^:]*\):\([0-9]*\):[[:space:]]*/\1: /' | sort > "$tmp"
if grep -vE '#\[ignore = "S[0-9]+' "$tmp" | grep -q .; then
  echo "untagged #[ignore] (must be #[ignore = \"S<n>: reason\"]):"; grep -vE '#\[ignore = "S[0-9]+' "$tmp"; exit 1
fi
if ! diff -u tests/IGNORED.txt "$tmp"; then echo "tests/IGNORED.txt is out of date; regenerate deliberately"; exit 1; fi
echo "ignore-drift: $(wc -l < "$tmp") ignored tests, all tagged, list matches"
