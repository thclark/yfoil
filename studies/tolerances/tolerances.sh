#!/usr/bin/env bash
# The tolerances study: measure the reference's own numerical noise floor.
#
# Every panel coordinate of the reference case is perturbed by +1 ULP (`math.nextafter`), the
# instrumented double-precision XFOIL reference is run on the base and on the perturbed panels,
# and the ES24.16 dumps of the two runs are compared line-aligned, stage by stage: lines,
# misaligned lines, max absolute difference, and the max, p99 and median relative difference
# (over values above 1e-8). Each run's iteration count and convergence are recorded beside it.
#
# This is a study, not a pipeline. Its findings were read, and the named constants in
# tests/common/utilities/tolerances.rs (TOL_PURE, TOL_LINALG, TOL_SOLVER, TOL_TRANSIENT) were
# set from them by hand (CLAUDE.md Rule 1); no test reads its output. Rerun it, and revisit
# those constants by hand, whenever the reference changes.
#
# Writes a run folder, studies/tolerances/runs/<UTC>/ (gitignored), whose publish/ holds the
# Markdown fragments the pages of docs/validation/tolerances/ include:
#   publish/tables/floor.md         the per-stage spread table
#   publish/tables/iterations.md    each run's iterations and convergence
#   publish/tables/run.md           which run this is
# It never writes under docs/: `cargo xtask publish tolerances` copies a run's publish/ there.
# The two runs work in target/noise-floor/{base,ulp}/. It builds the reference if it is missing.
#
# Usage: studies/tolerances/tolerances.sh [case-dir]   (default: tests/fixtures/xfoil/naca0012_n60_a2_re1e6)
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
CASE="${1:-$ROOT/tests/fixtures/xfoil/naca0012_n60_a2_re1e6}"
XF="$ROOT/target/xfoil-ref/instrumented/bin/xfoil"
[ -x "$XF" ] || "$ROOT/scripts/xfoil-build.sh"
W="$ROOT/target/noise-floor"; rm -rf "$W"; mkdir -p "$W/base" "$W/ulp"
python3 - "$CASE" "$W" <<'PY'
import math,sys
case,W=sys.argv[1],sys.argv[2]
src=open(f'{case}/panels.dat').read().split('\n')
base=[src[0]]; ulp=[src[0]]
for l in src[1:]:
    if not l.strip(): continue
    x,y=map(float,l.split())
    base.append(f" {x:.17e}  {y:.17e}"); ulp.append(f" {math.nextafter(x,math.inf):.17e}  {math.nextafter(y,math.inf):.17e}")
open(f'{W}/base/panels.dat','w').write('\n'.join(base)+'\n'); open(f'{W}/ulp/panels.dat','w').write('\n'.join(ulp)+'\n')
PY
for d in base ulp; do ( cd "$W/$d" && sed '/^CPWR\|^DUMP/d' "$CASE/xfoil.inp" | "$XF" > stdout.txt 2>&1 || true ); done
STAMP="$(date -u +%Y-%m-%dT%H-%M-%SZ)"
RUN="$ROOT/studies/tolerances/runs/$STAMP"; T="$RUN/publish/tables"; mkdir -p "$T"
DIRTY=""; [ -n "$(git -C "$ROOT" status --porcelain)" ] && DIRTY=", dirty"
echo "Results from the run of \`studies/tolerances/tolerances.sh\` at $STAMP on \`${CASE#"$ROOT"/}\` (git $(git -C "$ROOT" describe --tags --always)$DIRTY, host $(uname -m)-$(uname -s | tr '[:upper:]' '[:lower:]'))." > "$T/run.md"
python3 - "$W" "$T" <<'PY'
import re,sys
W,T=sys.argv[1],sys.argv[2]
sys.stdout=open(f'{T}/floor.md','w')
tok=re.compile(r'[-+]?(?:\d+\.\d*|\.\d+)(?:[EeDd][-+]?\d+)?')
def parse(p): return [[float(t.replace('D','E')) for t in tok.findall(l)] for l in open(p,errors='ignore')]
files=[('xfoil_panels.dat','geometry after ABCOPY'),('xfoil_inviscid.dat','inviscid GAM/QINV'),('xfoil_dij.dat','DIJ matrix'),
       ('mrchdu_input_1.dat','BL state after MRCHUE init'),('mrchdu_output_1.dat','BL state after MRCHDU #1'),
       ('blsolv_input.dat','SETBL VA/VB/VM/VDEL (calls 1..3)'),('blsolv_output.dat','BLSOLV VDEL out (calls 1..3)'),
       ('update_output_1.dat','UPDATE state/RLX/RMSBL (call 1)'),('viscal_iter.dat','per-iteration CL/CD/RMSBL/RLX (VISCAL call 1)'),
       ('viscal_final.dat','converged point (VISCAL call 1)'),('viscal_iters_all.dat','per-iteration record, all VISCAL calls'),
       ('viscal_points.dat','per-point CL/CD/CM/XTR (all VISCAL calls)')]
print(f"| file | lines | misaligned | max abs diff | max rel (abs>1e-8) | p99 rel | median rel | stage |")
print("|---|---|---|---|---|---|---|---|")
for f,desc in files:
    try: A=parse(f'{W}/base/{f}'); B=parse(f'{W}/ulp/{f}')
    except FileNotFoundError: print(f"| {f} | missing | | | | | | {desc} |"); continue
    n=min(len(A),len(B)); mis=sum(1 for i in range(n) if len(A[i])!=len(B[i])); absd=[]; rel=[]
    for i in range(n):
        if len(A[i])!=len(B[i]): continue
        for a,b in zip(A[i],B[i]):
            d=abs(a-b); absd.append(d)
            if max(abs(a),abs(b))>1e-8: rel.append(d/max(abs(a),abs(b)))
    rel.sort()
    if not rel: print(f"| {f} | {n} | {mis} | (no values) | | | | {desc} |"); continue
    print(f"| {f} | {n} | {mis} | {max(absd):.2e} | {rel[-1]:.2e} | {rel[int(0.99*(len(rel)-1))]:.2e} | {rel[len(rel)//2]:.2e} | {desc} |")
sys.stdout.close(); sys.stdout=open(f'{T}/iterations.md','w')
for d in ('base','ulp'):
    s=open(f'{W}/{d}/stdout.txt').read(); its=re.findall(r'^\s+(\d+)\s+rms:', s, re.M)
    print(f"\n{d}: iterations={its[-1] if its else '?'} converged={'VISCAL:  Convergence failed' not in s}")
sys.stdout.close()
PY
echo "run folder: $RUN"
echo "publish it with: cargo xtask publish tolerances --run $STAMP"
