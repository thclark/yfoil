# 7.7 TYPE 2 on a section that can reach negative CL drives itself onto the Kármán–Tsien pole — Replicated

**What TYPE 2 asks XFOIL to do.** TYPE 2 (`MATYP = 2`) is the fixed-lift type: it models an aircraft in
steady level flight, holding weight and chord fixed and letting the speed vary, so that the Mach number
follows the lift as `M∞ = M₁/√CL` (the three types are defined in
[the XFOIL reference](../../../user-guide/reference/xfoil-reference.md#the-oper-type-command-and-the-meaning-of-matyp-and-retyp)). It therefore does
not fix the Mach number. It asks for a state in which two conditions hold at once: the Mach number follows
the lift, `M∞ = M₁/√CL`, and the lift is whatever the aerofoil produces at that Mach, `CL = CL(α, M∞)`. SPECAL (`xoper.f`) solves this by Newton
iteration on a variable `CLM`, starting from `CLM = 1.0` on every call. That is the whole mechanism, so
the quality of the answer depends entirely on whether the fixed-point problem is well posed at the given α.

**Why the equation is well posed above zero lift and hostile below it.** The boundary is set by the
aerofoil's own compressibility correction. CLCALC forms Cp through the Kármán–Tsien denominator
`DEN = β + BFAC·Cp_inc`, which crosses zero when the local speed reaches `q/Q∞ = √(1 + 2β(1+β)/M∞²)`.
That is a **pole, not a loss of accuracy**: below it Cp is a smooth function of Mach, above it Cp has
changed sign (§4, CPCALC/CLCALC row).

At positive lift the coupling is **self-stabilising** — more lift lowers the Mach number, which weakens the
compressibility correction, which moves the denominator further from zero. At negative lift it is
**self-trapping**, in four steps:

1. `CL ≤ 0` puts the case outside the type's domain of definition — the speed at which a fixed weight is
   carried at negative lift is imaginary — so MRCL's floor `CLA = MAX(CLS, 1e-6)` stands in for it, and M∞
   is pushed to its *maximum* and clamped at 0.99.
2. At M∞ = 0.99 the pole sits at `q/Q∞ = 1.153`, which a 10 % section exceeds by α = ±1°.
3. Past the pole the denominator is negative, so **CL comes back positive at negative α**.
4. A positive CL is perfectly admissible to MRCL, so the fixed-point equation acquires a spurious root.
   Because `CL(α, M)` has a pole it takes every value in a neighbourhood of it, so a root always exists
   *adjacent to the pole*. SPECAL converges onto it in a handful of iterations and reports success.

**Method.** Measured on `naca64a010_n240_inviscid_polar30_type2_m03` (2026-09-17) by exact offline replay of
SPECAL's Newton. The replay is legitimate because on the inviscid path
`GAM(I) = cos α · GAMU(I,1) + sin α · GAMU(I,2)` is Mach-independent — Mach enters only through `BETA` and
`BFAC` inside CLCALC — so `CL(α, M)` can be re-evaluated for any M from the GAM already dumped at that α,
and the whole iteration replayed. The replay reproduces all 61 tracked points with zero relative
disagreement, so the figures below are XFOIL's own arithmetic.

**Every negative-α point is wrong**, not only the conspicuous ones, and most of them converge:

| α | XFOIL reports | truth (M→0) | Newton | min DEN |
|---:|---:|---:|---|---:|
| −30° | **+1.111** | −3.390 | converged, 6 it | −0.949 |
| −20° | **+1.040** | −2.318 | converged, 5 it | −0.0023 |
| −10° | **+0.327** | −1.176 | converged, 4 it | −0.0158 |
| −5° | **+0.128** | −0.590 | converged, 13 it | −0.361 |
| −1° | **−52.086** | −0.118 | failed, 20 it | −0.030 |

The Newton fails outright (`Minf convergence failed`) at only 6 of the 61 points (α = −11, −8, −3, −2, −1, 0).
The converged points are the worse case, because they carry no warning at all: XFOIL reports positive lift on
a symmetric aerofoil at −30° and prints nothing. α = 0 looks innocent only because CL = 0 is pinned by
symmetry whatever the Mach; its Newton failed just as hard, emitting 186 Mach-limit messages. α = +20° and
above is the same failure reached from the other side, once the suction peak crosses the M ≈ 0.3 threshold
of `q/Q∞ = 6.513`.

**The size of the anomaly tracks distance from the pole, not severity of the failure.** At α = −20° the
Newton converged to `M∞ = 0.294270`, where the minimum denominator is −0.0023 — essentially *on* the pole:

```
M = 0.294000   DEN = -0.000435   CL =   +59.89
M = 0.294270   DEN = -0.002316   CL =    +1.07   <- XFOIL's converged answer
M = 0.295000   DEN = -0.007411   CL =  -137.84
```

A change of 1e-3 in M∞ moves CL from +60 to −138, and XFOIL happens to land where it reads +1.04. By the same
token α = −1° reports −52.09 while α = −5° reports +0.128: the first is spectacular only because its
denominator is −0.030 rather than solidly negative, and the second is tamer for being *more* deeply inverted.
Both are equally broken. A polar plotted from this case therefore looks locally plausible exactly where it is
not, and a departure metric built on the size of the anomaly will rank these points in the wrong order.

**These points are ill-conditioned pointwise, not by accumulation** (`docs/development/conventions/terminology.md`). They reproduce bitwise,
and forward (−30→+30) and reverse (+30→−30) sweeps agree at all 61 points, because SPECAL resets `CLM = 1.0`
at every call and GAMU is α-independent: nothing carries between points on the inviscid path. Their wide
1-ULP twin envelope is pointwise ill-conditioning — each sits next to a pole, so a 1-ULP perturbation
genuinely does move it by O(1). The case is consequently a useful control for sensitivity studies: a wide
envelope that is provably non-cumulative. It says nothing either way about accumulation on the viscous
cases, where BL state really is carried forward.

**yFoil.** Replicated: `set_mach_re_from_cl` (`src/solver/setbl.rs:40-82`) carries both MRCL guards and
`clcalc.rs:57-69` the unguarded denominator, so yFoil converges onto the same spurious roots — and
then flags them, rather than reporting the numbers as if they were sound. This case is the worked
example behind [Solution validity](../../../user-guide/reference/solution-validity.md), which catalogues every reason yFoil withholds
a result and what each one means. **Test:** `tests/known_issues/kt_pole_7_7.rs` runs the same section at N = 60 at
−23°, 0° and 19° (`naca64a010_n60_inviscid_type2_m03_am23_a0_a19`: the three angles take every branch, and raise
every event, of the −30…30° sweep, the points being independent) and checks yFoil's validity record against
XFOIL's own events point by point. The measurements above were made at N = 240 (`naca64a010_n240_inviscid_polar30_type2_m03`,
defined in `cases.toml` and regenerated by `cargo xtask fixtures --case` with that name; untracked). The
negative leg of such a sweep must not be read as a polar, and it is excluded from any accuracy claim.

## Runs that show it

- [NACA 64A010 inviscid, TYPE 2 at M 0.3](../runs/naca-64a010-inviscid-type-2-at-m-0-3.md)
- [NACA 64A010, TYPE 2 at M 0.3](../runs/naca-64a010-type-2-at-m-0-3.md)
