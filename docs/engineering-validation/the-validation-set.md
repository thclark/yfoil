# The validation set

The sections and flow conditions the rest of this tab assesses, and why each was chosen. Each
section is public, mathematically defined and reproducible from its designation; its definition
is in [Aerofoil families](../using-yfoil/aerofoil-families/index.md).

## The sections

One section per generator family (two from the 4-digit series), each chosen for something the
others do not exercise.

| Section | Family | Why | Generate |
|---|---|---|---|
| NACA 0012 | [4-digit](../using-yfoil/aerofoil-families/naca-4-digit-and-4-digit-modified.md) | Symmetric baseline and CI reference case: the XFOIL-independent invariants (C_L = C_M = 0 at α = 0, mirror symmetry at −α) and every solver stage's first fixture. ([jacobs1933](../using-yfoil/bibliography.md#jacobs1933)) | `yfoil geometry naca 0012` |
| NACA 4412 | [4-digit](../using-yfoil/aerofoil-families/naca-4-digit-and-4-digit-modified.md) | The cambered 4-digit member, a real wing section, and the ±15° polar case that runs past the limits of convergence. ([jacobs1933](../using-yfoil/bibliography.md#jacobs1933), [lednicer](../using-yfoil/bibliography.md#lednicer)) | `yfoil geometry naca 4412` |
| NACA 0012-34 | [4-digit modified](../using-yfoil/aerofoil-families/naca-4-digit-and-4-digit-modified.md#naca-4-digit-modified) | The 0012 thickness with a sharper nose (I = 3) and its maximum thickness aft (M = 0.4): leading-edge radius separated from thickness. ([stack1934](../using-yfoil/bibliography.md#stack1934)) | `yfoil geometry naca 0012-34` |
| NACA 23018 | [5-digit](../using-yfoil/aerofoil-families/naca-5-digit-including-reflex.md) | A real aircraft root section (Beech 90 King Air) with the camber concentrated at the nose: transition and separation are exercised at the nose rather than mid-chord. ([jacobs1935](../using-yfoil/bibliography.md#jacobs1935), [lednicer](../using-yfoil/bibliography.md#lednicer)) | `yfoil geometry naca 23018` |
| NACA 16-212 | [16-series](../using-yfoil/aerofoil-families/naca-16-series.md) | Designed for high subsonic speed and tested at Mach 0.3–0.8: the compressible (Kármán–Tsien) case. ([stack1943](../using-yfoil/bibliography.md#stack1943), [lindsey1948](../using-yfoil/bibliography.md#lindsey1948)) | `yfoil geometry naca 16-212` |
| NACA 63-415 | [6-series](../using-yfoil/aerofoil-families/naca-6-and-6a-series.md) | The classic wind-turbine section, with published XFOIL and Navier–Stokes comparisons; a closed trailing edge; the family without a closed-form generator. ([abbott1945](../using-yfoil/bibliography.md#abbott1945), [bertagnolio2001](../using-yfoil/bibliography.md#bertagnolio2001)) | `yfoil geometry naca 63-415` |
| NACA 64A010 | [6A-series](../using-yfoil/aerofoil-families/naca-6-and-6a-series.md#naca-6a-series) | The standard transonic test section: symmetric, straight-sided aft, sharp trailing edge with a finite angle. ([loftin1948](../using-yfoil/bibliography.md#loftin1948)) | `yfoil geometry naca 64A010` |
| Kármán–Trefftz (xc = −0.1, yc = 0.05, τ = 10°) | [analytic](../using-yfoil/aerofoil-families/karman-trefftz-sections.md) | Exact potential-flow solution: the panel method checked against a closed-form answer, and a sharp trailing edge on a smooth, non-NACA shape. ([karman1918](../using-yfoil/bibliography.md#karman1918), [joukowski1910](../using-yfoil/bibliography.md#joukowski1910)) | `yfoil geometry karman-trefftz --x-centre -0.1 --y-centre 0.05 --te-angle 10` |

## The flow conditions

The comparison with XFOIL runs these sections through the regimes people use them in, under
normal conditions: Reynolds numbers from 1e4 to 1e8, Ncrit from 4 (a disturbed stream) to 11 (a
quiet one), Mach numbers from 0 to 0.6, transition forced by trips, a fixed-CL sweep, and an
inviscid sweep of the analytic section. The polars are swept from 0° to ±30° by XFOIL's own polar
procedure, at 160 panels and up to 100 iterations per point, which represents the finer
panellings wherever the flow is attached
([known issues §7.13](../fidelity-to-xfoil/xfoils-known-issues/the-issue-register.md)). Every case,
with its conditions, is listed and plotted on
[agreement with XFOIL under normal conditions](agreement-with-xfoil-under-normal-conditions.md);
the cases are `group = "series"` in `xtask/fixtures-config/cases.toml`.

Runs chosen to show one of XFOIL's known issues — the dependence on the panel count, the
reference hanging on a non-finite value — are deliberately *not* in this set: they are not normal
conditions, and they are shown under Fidelity to XFOIL as
[runs that show the issues](../fidelity-to-xfoil/xfoils-known-issues/runs-that-show-the-issues.md).
