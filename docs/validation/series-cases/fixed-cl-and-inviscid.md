---
icon: lucide/crosshair
---

# Fixed CL and inviscid

The two other ways of driving OPER: a fixed-CL sweep (one `CL` command per point, so SPECCL and MRCL run at every point), and an inviscid sweep of the Kármán–Trefftz section (one SPECAL per degree, no boundary layer).

--8<-- "validation/series-cases/tables/cases-fixed-cl-and-inviscid.md"

## NACA 4412, fixed CL

`naca4412_n160_clsweep_re1e6_iter100` — NACA 4412, fixed CL from -0.5 to 1.5 by 0.1: SPECCL and MRCL at every point.

<ironlab-figure src="figures/polars-naca4412_n160_clsweep_re1e6_iter100.fig" format="fig" name="polars-naca4412_n160_clsweep_re1e6_iter100" alt="naca4412_n160_clsweep_re1e6_iter100: CL, CD, the last RMSBL, iterations, transition location and position in its station interval against alpha, XFOIL and yFoil">naca4412_n160_clsweep_re1e6_iter100: CL, CD, the last RMSBL, iterations, transition location and position in its station interval against alpha, XFOIL and yFoil</ironlab-figure>

!!! note "Stub — to be written"

    What this run shows, and where and why yFoil and XFOIL part.

## Kármán–Trefftz, inviscid

`kt_n160_inviscid_polar30` — Kármán–Trefftz, the sharp analytic section, inviscid, one SPECAL per°ree.

<ironlab-figure src="figures/polars-kt_n160_inviscid_polar30.fig" format="fig" name="polars-kt_n160_inviscid_polar30" alt="kt_n160_inviscid_polar30: CL, CD, the last RMSBL, iterations, transition location and position in its station interval against alpha, XFOIL and yFoil">kt_n160_inviscid_polar30: CL, CD, the last RMSBL, iterations, transition location and position in its station interval against alpha, XFOIL and yFoil</ironlab-figure>

!!! note "Stub — to be written"

    What this run shows, and where and why yFoil and XFOIL part.
