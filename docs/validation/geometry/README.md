# Geometry validation

Two generators, two questions.

**What matters for solver equivalence** is that YFoil's panels survive the hand-off into XFOIL
bitwise (CLAUDE.md Rule 4): YFoil writes `panels.dat` at 17 significant digits, XFOIL `LOAD`s it,
and the fixture pipeline asserts XFOIL's `X/Y` after `ABCOPY` are identical to what was written
(`cargo xtask fixtures`, every case, every run). That has held on all cases since 2026-09-03.

**Stage G** is the optional question: can YFoil reproduce what XFOIL itself produces from a bare
`NACA dddd` command? That needs XFOIL's own generator (`naca.f`, thickness applied vertically,
NSIDE = 123 points per side, TE-bunched spacing with AN = 1.5) and `PANGEN`. Both are translated
line for line (`naca_4digit_xfoil`, `naca_5digit_xfoil`, `repanel_xfoil`; CLI `yfoil geom naca
--naca-model xfoil`).

## Gate (2026-09-04)

Reference: instrumented XFOIL 6.99 (double precision), `NACA dddd / PPAR / N n`, PANGEN dumping
the buffer airfoil `XB/YB/SB`, the paneled `X/Y/S`, `SBLE`, `SLE` and `SHARP` at `ES24.16`
(`tests/fixtures/xfoil/pangen_*/xfoil_pangen.dat`). Test: `tests/xfoil_pangen_tests.rs`.

| case | buffer points | buffer values bitwise identical | PANGEN nodes | nodes bitwise identical | worst node \|diff\| |
|---|---|---|---|---|---|
| NACA 0012, NPAN 160 | 245 | 735 / 735 | 160 | 106 / 480 | 6.7e-16 |
| NACA 0012, NPAN 81 | 245 | 735 / 735 | 81 | 42 / 243 | 8.9e-16 |
| NACA 4412, NPAN 160 | 245 | 735 / 735 | 160 | 38 / 480 | 1.6e-15 |
| NACA 4412, NPAN 81 | 245 | 735 / 735 | 81 | 48 / 243 | 8.3e-16 |
| NACA 23012, NPAN 160 | 245 | 735 / 735 | 160 | 96 / 480 | 8.9e-16 |

The buffer airfoils are bit-identical (same `pow`/`sqrt` on the same host). The PANGEN nodes
differ at the last ULP or two — the spline evaluations in the node-placement Newton iteration
are associated slightly differently by the two compilers — and are gated at `TOL_PURE` (1e-12).
`SBLE`, `SLE` and the `SHARP` flag match.

## What the old table measured

The previous version of this page reported RMS errors of 4.9e-4 (0012) to 1.6e-3 (4412). Two
things were stacked in those numbers:

1. **A 7-digit reference.** The comparison used XFOIL's `PSAVE` output, which is
   `FORMAT(1X,G15.7,G15.7)`: seven significant figures. Nothing below ~1e-7 could ever be shown.
2. **Two different NACA definitions.** XFOIL's `NACA4` applies the thickness distribution
   vertically (`YB = YC ± YT`, `naca.f:62`). The NACA definition applies it perpendicular to the
   camber line, which is what YFoil's default generator (`naca_4digit`) does — a deliberate,
   documented divergence (CLAUDE.md Rule 2). For a symmetric section the two agree (the camber
   slope is zero), which is why 0012 was "better" than 4412 by 3×: the 4412 figure was the
   thickness-application difference, not a paneling error.

Neither generator is "wrong"; they answer different questions. The default stays the exact
NACA definition; `--naca-model xfoil` reproduces XFOIL's.

## Reproducing

```bash
cargo xtask fixtures --case pangen_naca0012_n160 --case pangen_naca4412_n81   # regenerate the dumps
cargo test --test xfoil_pangen_tests -- --nocapture                            # the gate, with counts
yfoil geom naca 4412 -n 160 --naca-model xfoil -o naca4412_xfoil.json          # XFOIL's panels from YFoil
```
