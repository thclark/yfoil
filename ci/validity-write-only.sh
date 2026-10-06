#!/usr/bin/env bash
# The validity record is observational (CLAUDE.md Rule 2): `src/solver/` and `src/bl/` write it
# and must never read it back, so that no numeric result can depend on a validity classification.
# Every mention of `.validity` in those trees must be a write — an assignment to one of its
# fields — or the definition itself. A read (a condition, a comparison, a call to one of the
# classifier methods) fails here rather than silently becoming a control decision.
set -euo pipefail
cd "$(dirname "$0")/.."

# the record's own definition and the state field that holds it are not solver logic
EXCLUDE='src/solver/validity.rs|src/solver/blstate.rs|src/solver/analysis.rs'

hits=$(grep -rn '\.validity' src/solver src/bl --include='*.rs' | grep -vE "^($EXCLUDE):" || true)

fail=0
while IFS= read -r line; do
  [ -z "$line" ] && continue
  code=${line#*:}; code=${code#*:}
  # a write is `<expr>.validity.<field> = <expr>;` with no comparison on the left of the `=`
  if echo "$code" | grep -qE '^\s*[A-Za-z_][A-Za-z0-9_]*\.validity\.[a-z_]+ = [^=]'; then
    continue
  fi
  echo "VALIDITY-READ: $line"
  fail=1
done <<< "$hits"

if [ $fail = 0 ]; then
  echo "validity-write-only: clean ($(echo "$hits" | grep -c . || true) writes in src/solver, src/bl)"
fi
exit $fail
