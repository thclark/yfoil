# 7.13 Whether a sweep converges, and what it converges to, can depend on the panel count — Replicated

The same OPER script on the same section, panelled at N = 60, 160 and 240 nodes (yFoil's cosine
sampling; ITER 20, 100 and 200, the settings each was studied at), can differ qualitatively, not by
a discretisation error that shrinks with N. Measured 2026-10-08 (`cargo xtask fixtures`, reference
runs; ITER 200 at N = 160 as a control, so the difference is the panelling's and not the iteration
limit's):

| Case (0 → +30°, `INIT`, 0 → −30° by 1°) | N = 60, ITER 20 | N = 160, ITER 100 | N = 160, ITER 200 | N = 240, ITER 200 |
|---|---|---|---|---|
| NACA 16-212, Re 1e6, M 0.7 | first unconverged point 13° | first unconverged 3°; non-finite at −18° and hangs (§4) | first unconverged 3°; runs to its end | first unconverged 23°; CL_max 1.257 at 26° |
| NACA 64A010, Re 1e6, M 0, TYPE 2 | 5 of 14 points converged | 11 of 13 converged, CL to 0.95; non-finite at 12° and hangs | 11 of 15 converged; non-finite at 14° and hangs | **no point converged** (10 recorded, every one at the limit) |
| NACA 4412, inviscid, M 0.7 | largest CL 111 (at −15°) | largest CL 18 (16°) | — | largest CL 73 (18°) |
| NACA 64A010, inviscid, M 0.3, TYPE 2 | largest CL 36 (−1°) | largest CL 18 (−1°) | — | largest CL 28 (−3°) |

The inviscid rows are the Kármán–Tsien pole (§4, §7.7): past it the corrected Cp, and so CL, depend
on which panel's speed crosses the threshold first and by how much, so their size is the
panelling's choice. The viscous rows are convergence itself: N = 240 does not converge where N = 160
does, and the reverse at low alpha on the 16-212. Neither panel count is a reference truth.

Where the sections are attached and well conditioned, N = 160 and N = 240 agree to a few per cent in
CL and CD and on CL_max and its angle (NACA 0012 and 4412 at Re 1e6, the 4412 at Re 3e6 and M 0.6,
the Kármán–Trefftz section, the fixed-CL sweep); past stall the branch a sweep follows can differ
between them, which is expected. That is why the series cases run at N = 160, ITER 100, and these
four runs are known-issues cases at all three panel counts (`group = "known-issues"`,
`docs/validation/known-issues/`). yFoil, given the same panels, is to follow XFOIL at each.

## Runs that show it

- [NACA 16-212 at M 0.7](../runs/naca-16-212-at-m-0-7.md)
- [NACA 64A010, TYPE 2](../runs/naca-64a010-type-2.md)
- [NACA 4412 inviscid at M 0.7](../runs/naca-4412-inviscid-at-m-0-7.md)
- [NACA 64A010 inviscid, TYPE 2 at M 0.3](../runs/naca-64a010-inviscid-type-2-at-m-0-3.md)
