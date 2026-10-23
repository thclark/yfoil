---
icon: lucide/hourglass
---

# 9. Pathological cases

Cases whose reference runs are pathologically slow or hang, kept as the `pathological` group of
`xtask/fixtures-config/cases.toml` (untracked, run only by name or `--group pathological`, never by
`--big`, never in a routine regeneration). Each was a copy of one of the removed twins study's N = 240,
ITER 200 cases. They are recorded here with what the run showed and what physics it captures, so
the insight is kept even though the run is not part of the routine sets; the last column says where
the same behaviour is tested and shown in routine time.

| | Case | What makes it pathological (observed 2026-10-07) | Physical insight / scenario | Where the behaviour is tested and shown |
|---|---|---|---|---|
| 9.1 | `pathological_naca4412_n240_polar30_re1e6_iter200` (NACA 4412, Re 1e6, ±30°) | Past about 30° on the upward leg CL turns non-finite while the boundary-layer state stays finite, and XFOIL hangs formatting the plot label (§4); the downward leg never runs. Twins about 23 min each. | The cambered section through stall: the reference cannot complete a post-stall sweep at this panelling and iteration limit; its twins leave the reference's branch at the stall edge. | Tested: the case-cover member `naca4412_n160_polar_down16_re1e6_iter100` (whole run to call 30, every iteration of call 31) and the tracked N = 60 sweep's steps. Shown: series `naca4412_n160_polar30_re1e6_iter100` (at N = 160 the reference does not hang); the hang itself, reference only: known-issues `naca4412_n240_polar30_re1e6_iter200`. |
| 9.2 | `pathological_naca4412_n240_polar30_re3e6_iter200` (NACA 4412, Re 3e6, ±30°) | Hangs in the plot label like the Re 1e6 sweep. Twins 80+ min each; one twin goes non-finite at a point the reference converges (CL envelope 0.78). | MRCHUE's iterate fallback on the first march at high alpha: at Re 3e6 the fresh march stalls, and the post-stall solution is ill-conditioned. | Tested: case-cover member `naca4412_n60_a18_re3e6_iter40` (every iteration of call 1). Shown: series `naca4412_n160_polar30_re3e6_iter100`; the hang itself, reference only: known-issues `naca4412_n240_polar30_re3e6_iter200`. |
| 9.3 | `pathological_naca0012_n240_polar30_re1e6_type3_iter200` (NACA 0012, Re 1e6, TYPE 3, ±30°) | Twins up to about 15 min each (reference about 4 min); the twins leave the reference's branch from 17–20° and at −29…−30° (6 completion mismatches, CL envelope 8e-2). | TYPE 3 scales Re with 1/CL: near CL = 0 the Reynolds number runs to MRCL's limit (`MRCL_RE_LIMIT` fires at every call), and in stall the solution is strongly ill-conditioned. | Tested: case-cover member `naca0012_n60_a2_re1e6_type3` (whole run and a step-cover step). Not shown by a series case (removed 2026-10-08 as pathological). |
| 9.4 | `pathological_naca23012_n240_polar30_re1e6_xtr022_0001_iter200` (NACA 23012, Re 1e6, trips 0.22 / 0.001, ±30°) | Non-finite from the 2° point; every station Newton then fails, about 7 s per SETBL call, and the full sweep runs for hours. | **The lower trip at x/c = 0.001 lies upstream of the stagnation point at positive alpha** (XIFSET's "past trip"): TRCHEK2 is asked for transition upstream of a laminar station, produces a NaN near the leading edge, and the march continues on the fallbacks. Reached from the 1° solution, at N ≥ 100, it takes TRCHEK2's `IF(AX .LE. 0.0)` exit (`xblsys.f:390`), once annotated unreachable (§8). | Tested: case-cover member `naca23012_n60_a4_re1e6_xtr022_0001`, and the non-finite probe `naca23012_n100_a1_a2_re1e6_xtr022_0001` (the `xblsys.f:390` exit: a step-cover step and every iteration of call 2). Shown: known-issues `naca23012_n60_polar30_re1e6_xtr022_0001` (non-finite from 1°, the trip-upstream event at every station of every call after it). |

The removed `…_iter200_stall` copy repeated the first row's conditions exactly and is not kept.

**Slow only through its twins: run time as a symptom of ill-conditioning.** Three of the four cases
are not slow in their own reference run; their seeded 1-ULP twins are — the same script, every panel
coordinate moved by −1, 0 or +1 ULP. Measured 2026-10-07 (wall clock; the reference ran alone, its
five twins together on 8 cores, so some of a ratio can be load, but not a factor of ten or more):

| Case | Reference | Its five twins | Ratio |
|---|---|---|---|
| NACA 4412, Re 1e6 | 132 s, of which 120 s is the watchdog's wait at the plot-label hang (about 12 s of computation) | 1302–1423 s each | about 10× on wall time, about 100× on computation |
| NACA 4412, Re 3e6 | 160 s, of which 120 s is the wait at the hang (about 40 s of computation) | 4855–5134 s each | about 30× on wall time, over 100× on computation |
| NACA 0012, TYPE 3 | 252 s | 236–874 s (three twins as fast as the reference, two 3.5× slower) | 1–3.5× |
| NACA 23012, trips | hours in its own run | — | — |

A one-ULP change to the geometry changes the solution time by up to two orders of magnitude: past
stall the twins leave the reference's branch (the completion mismatches above) and take another route
through the sweep. Which route, and why it costs so much more — more points iterated to their limit,
marches on the fallbacks, or not reaching the non-finite CL that ends the reference at the hang — is
for the later study to establish (§8). Solution time is therefore itself a measure of how ill-conditioned a
run is, and an interesting case in itself for the later study: the same physics, identical to the last
bit but one, at a hundred times the cost.

## Runs that show it

- [NACA 23012 with a trip upstream of stagnation](../runs/naca-23012-with-a-trip-upstream-of-stagnation.md)
