---
icon: lucide/trending-up
---

# Polars

``` sh
yfoil polar <GEOMETRY.json> [OPTIONS]
```

``` sh
yfoil polar naca4412.json \
    --alpha-min -5 --alpha-max 15 --alpha-step 0.5 \
    --reynolds 1e6 --max-iterations 20 \
    --label "NACA 4412" --distributions -o naca4412_polar.json
```

## How the sweep runs

A polar is a **state machine**, not N independent solves: the boundary layer
converged at one angle is the initial condition for the next. yFoil follows
XFOIL's procedure exactly:

1. Solve 0°, then sweep **0° → α<sub>max</sub>**.
2. Reinitialise the boundary layer and re-solve 0°.
3. Sweep **0° → α<sub>min</sub>**.
4. Stitch the two legs into one ascending polar.

The re-solved 0° seeds the downward leg and is not itself a polar point. There is
no configurable start angle — both legs start from 0°, as XFOIL does.

## Options

| Option | Default | Meaning |
| --- | --- | --- |
| `--alpha-min` | -5 | Lowest angle of attack, degrees |
| `--alpha-max` | 15 | Highest angle of attack, degrees |
| `--alpha-step` | 0.5 | Step size, degrees |
| `-r`, `--reynolds` | 1e6 | Reynolds number |
| `-m`, `--mach` | 0 | Mach number |
| `--ncrit` | 9 | Critical amplification factor |
| `--max-iterations` | 20 | Maximum viscous Newton iterations per point |
| `--label` | file stem | Legend entry for `yfoil plot polar` |
| `--distributions` | off | Embed the full analysis record at every point |
| `--include-lagged-closures` | off | Also embed the lagged closure arrays |
| `--json` | off | Emit JSON on stdout |
| `-o`, `--output` | — | Output file |

`--distributions` makes the polar file self-contained: it can then be fed to
`yfoil plot foil` to draw boundary-layer distributions at selected angles.

## The XFOIL equivalent

``` text
OPER
VISC 1000000
ITER 20
PACC
polar_pos.txt

ALFA 0
ASEQ 1 15 1
PACC

INIT
PACC
polar_neg.txt

ALFA 0
ASEQ -1 -15 -1
```
