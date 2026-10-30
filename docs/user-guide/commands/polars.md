---
icon: lucide/trending-up
---

# Polars

``` sh
yfoil polar <GEOMETRY.json> [OPTIONS]
```

``` sh
yfoil geometry naca 4412 -n 160 -o naca4412.json
yfoil polar naca4412.json \
    --alpha-min -5 --alpha-max 15 --alpha-step 0.5 \
    --reynolds 1e6 --max-iterations 20 \
    --label "NACA 4412" --distributions -o naca4412_polar.json
```

Drawn with [`yfoil plot polar`](plotting.md#polar) and [ironlab](https://ironlab.org),
this sweep produces the polar below.

<ironlab-figure src="../figures/naca4412_polar.fig" format="fig" name="naca4412_polar" alt="NACA 4412 polar at Re 1e6, alpha -5 to 15 degrees">NACA 4412 polar at Re 1e6, alpha -5 to 15 degrees</ironlab-figure>

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
| `--allow-invalid` | off | Show the numbers of invalid points of any class (out-of-domain, unconverged, ill-conditioned); their status is unchanged |
| `--conditioning` | `on` | `on`: check each point with its five 1-ULP twins. `off`: run no twins, as XFOIL. `study`: `on`, and also write every twin's record under `study` |
| `--halt-on` | `xfoil` | Halt a leg after N consecutive points meeting a criterion: `CRITERION[=N],...` ([below](#halting-a-leg)) |
| `--spread-warning-threshold` | 1e-8 | Spread over the twins beyond which a valid point carries a warning that it may be ill-conditioned |
| `--spread-error-threshold` | 1e-4 | Spread over the twins beyond which a point is ill-conditioned: `Invalid`, numbers withheld |
| `--json` | off | Emit JSON on stdout |
| `-o`, `--output` | — | Output file |

`--distributions` makes the polar file self-contained: it can then be fed to
`yfoil plot foil` to draw boundary-layer distributions at selected angles.

Each point is checked with its [twins](../reference/solution-validity.md#twins) unless
`--conditioning off` is given. A twin of a polar point starts from the same boundary layer the point starts
from, its predecessor's, so the check is of that one solve and does not depend on how the sweep got
there. Points beyond the warning threshold only are listed `Valid` with a warning that they may
be ill-conditioned. Ill-conditioned points are listed as `Invalid (IllConditioned)` with their numbers withheld,
counted in the summary, and left out of \(C_{L,max}\), \((L/D)_{max}\) and \(C_{D0}\); the
sweep itself runs exactly as XFOIL's, since the twins never steer it. `--conditioning study` keeps
every twin's record ([below](#keeping-the-twins-studies)).

## Halting a leg

XFOIL halts a sequence after four consecutive points whose viscous iteration did not converge
(`NSEQEX`), and yFoil does the same by default. `--halt-on` takes a list of criteria, each with
the number of consecutive points that halts a leg:

| Criterion | Counts a point that is | Default N |
| --- | --- | --- |
| `xfoil` | viscous and unconverged (XFOIL's own rule) | 4 |
| `out-of-domain` | invalid for a Class A reason | 1 |
| `unconverged` | invalid for a Class B reason | 1 |
| `ill-conditioned` | invalid for a Class C reason (needs the twins) | 1 |
| `none` | — never halt | — |

Each criterion keeps its own count, reset by a point that does not meet it, so
`--halt-on xfoil,ill-conditioned` keeps XFOIL's rule and also stops a leg at its first
ill-conditioned point — typically past stall, before the solutions diverge. `xfoil` is narrower
than `unconverged`: XFOIL does not count the inviscid Newtons, nor a converged point that passed
through a NaN, so only `xfoil` reproduces XFOIL's sweep exactly. The alphas a halted leg never
reaches are reported `NotAttempted`, and the polar's `halts` names the criterion and where it
stopped.

A NACA 4412 at Re 3e6, swept down to −25°, has its first ill-conditioned point at −18°:

``` sh
yfoil geometry naca 4412 -n 160 --method cosine -o naca4412.json
yfoil polar naca4412.json --alpha-min -25 --alpha-max 0 --alpha-step 1 --reynolds 3e6 \
    --max-iterations 100 --halt-on xfoil,ill-conditioned -o naca4412_halt.json
```

``` text
  -25.00          -          -          -        -        -        -     -          -  NotAttempted (SequenceHalted)
     ...
  -19.00          -          -          -        -        -        -     -          -  NotAttempted (SequenceHalted)
  -18.00          -          -          -        -        -        -     -          -  Invalid (IllConditioned); warning: solution converged but 1 twin did not
  -17.00   -1.34185   0.025649   -0.10173   -52.32    1.000    0.009     6    2.26e-5  Valid
     ...
Summary:
  ...
  Valid: 18/26 points
  Invalid: 1 point(s) — solved, but the numbers are not to be believed (numbers withheld; --allow-invalid reports them)
    of which ill-conditioned: 1 — a 1-ULP jitter of the panel coordinates moves them beyond 1e-4
  Not attempted: 7 point(s) — the sweep halted before reaching them

Note: a leg halted after 1 consecutive `ill-conditioned` point(s), at alpha = -18.00° (--halt-on)
```

The upward leg is untouched. In the JSON, `completed` is `false` and `halts` records the stop:

``` json
"halts": [
  { "criterion": "ill-conditioned", "after": 1, "alpha_deg": -18.0 }
]
```

## Keeping the twins: studies

`--conditioning study` keeps every twin's record, in a `study` array beside `results`, so that
everything the check computed can be traced. Run the sweep above with it, and the −18° point
(`"id": "bafgwwmm0sgs"`) has five twins in `study`. Three reach \(C_L = -0.778\), one
\(C_L = -0.662\), and one runs out of iterations at \(C_L = -1.241\); that spread of 0.37 in
\(C_L\) is why the point is ill-conditioned. The twin that did not converge is itself invalid, and
still carries its numbers:

``` json
{
  "alpha_deg": -18.0,
  "id": "92r0whuhgf9w",
  "initialised_from": "y76br8dein8n",
  "status": "Invalid",
  "reasons": ["ViscousNotConverged"],
  "classes": ["unconverged"],
  "role": "twin",
  "twin_of": "bafgwwmm0sgs",
  "seed": 2,
  "diagnostics": { "converged": false, "iterations": 105, "residual": 0.8830497281847924, ... },
  "values": { "cl": -1.2412261265650741, "cm": -0.07961596507156694, "cd": 0.051776523825930273, ... }
}
```

Each twin has an `id` of its own and cites, as `initialised_from`, the point's own predecessor: it
is the same solve from the same state, on jittered panels. A study's records always carry their
numbers, whatever their status and whatever `--allow-invalid` says; `results` is unchanged by a
study. With `--distributions`, each embedded record also carries its twins' full records
(geometry, surface and boundary layer) under its own `study`.

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
