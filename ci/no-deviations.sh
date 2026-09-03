#!/usr/bin/env bash
# Fail if src/ contains the vocabulary of temporary algorithm changes (CLAUDE.md Rule 2),
# except for lines covered by ci/deviations-allowlist.txt.
set -euo pipefail
cd "$(dirname "$0")/.."
PATTERNS='disable_vm|disable_dij|disable_coupling|TEMPORARILY DISABLED|coupling_scale|causes divergence|tighter than XFOIL|simplified for stability|max_iter_debug|for stability'
hits=$(grep -rnE "$PATTERNS" src --include='*.rs' || true)
fail=0
while IFS= read -r line; do
  [ -z "$line" ] && continue
  file=${line%%:*}; rest=${line#*:}; rest=${rest#*:}
  ok=0
  while IFS= read -r allow; do
    case "$allow" in ''|'#'*) continue;; esac
    a=${allow%%#*}; a=$(echo "$a" | sed 's/[[:space:]]*$//')
    afile=${a%%:*}; apat=${a#*:}
    if [ "$afile" = "$file" ] && echo "$rest" | grep -qF -- "$apat"; then ok=1; break; fi
  done < ci/deviations-allowlist.txt
  if [ $ok = 0 ]; then echo "DEVIATION: $line"; fail=1; fi
done <<< "$hits"
[ $fail = 0 ] && echo "no-deviations: clean (allowlist has $(grep -cvE '^\s*(#|$)' ci/deviations-allowlist.txt) tolerated entries)"
exit $fail
