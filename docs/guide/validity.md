---
icon: lucide/shield-check
---

# Solution validity

XFOIL will hand back numbers for an operating point it has not actually solved. Sometimes it says
something on the console first; sometimes it says nothing at all. Either way the numbers arrive
looking exactly like good ones, and a polar plotted through them looks like a polar.

yFoil reproduces XFOIL's arithmetic exactly, including the cases above — that is the point of the
project ([Rule 2](../xfoil-known-issues.md)) — so it computes the same wrong numbers. What it does
differently is **say so**, and withhold them unless you ask for them anyway.

## The contract

Every operating point you ask for comes back as a record. A point is never silently missing from a
result, and the record always carries a `status`:

```json
{
  "alpha_deg": -1.0,
  "id": "5273jwcqarap",
  "initialised_from": "y8v0z1ll6a2m",
  "status": "Invalid",
  "reasons": ["KarmanTsienOutOfDomain", "MachClNewtonExhausted"],
  "diagnostics": { "karman_tsien_margin_forces": -0.0296, "mach_limited": true },
  "values": null
}
```

`values` — the lift, drag and moment — is present only when the status is `Valid`. Pass
`--allow-invalid` to `yfoil analyse` or `yfoil polar` and it is filled in for invalid points too;
that flag changes **only** whether the numbers are present, never the status or the reasons, so a
result carrying numbers can still be recognised as one that should not be trusted.

For `yfoil analyse` the same gate covers `surface` (the q and Cp distributions) and
`boundary_layer`, because those are results too — on an out-of-domain point the Cp field is
precisely what went wrong. `conditions` and `geometry` are *inputs*, so they are always present.

Both commands **exit zero** whenever the solver ran to completion. An invalid point is a result,
not a tool failure, and the status already carries it; a non-zero exit is reserved for the tool
failing to do what was asked, such as unreadable geometry.

Summary statistics (`cl_max`, `ldratio_max`, `cd0`) and every plotted curve are computed from the
valid points alone, so an artefact cannot become the reported CL~max~.

## Status codes

| `status` | Meaning |
|---|---|
| `Valid` | Solved, every formula was evaluated inside its domain, and every iteration met its tolerance. `values` present. |
| `Invalid` | Solved, but at least one `reason` applies, so the numbers are not to be believed. `values` withheld unless `--allow-invalid`. |
| `NotAttempted` | Never solved. The sweep halted before reaching this alpha (XFOIL's NSEQEX rule ends a sequence after four consecutive non-converged points). It is reported because you asked for it. |

## How failures are classified

Three distinct things can be wrong with a solved point, and they differ in something more important
than severity: **what it would take to know**. Keeping them apart is what stops the scheme from
either over-claiming or quietly ignoring a real problem.

| | Class | What it is | Provable from | Shape |
|---|---|---|---|---|
| **A** | Domain violation | A formula was evaluated outside its analytic domain, so its result is wrong rather than imprecise | a single run | binary, local, with a continuous margin |
| **B** | Iteration exhaustion | An iteration hit its cap without meeting its tolerance, so the value returned is not a solution to the problem posed | a single run | binary, global to the point |
| **C** | Conditioning | The point is converged and in domain, but sits where a 1-ULP change of input moves it by O(1) | **two or more runs** | continuous; no threshold |

**A and B are what this page catalogues.** Both are exact facts about the computation as it was
performed, so they can be recorded as it runs, and both are grounds for withholding the numbers.
They are detected by observation only: the solver records the evidence and never reads it back, so
no numeric result can depend on a validity classification
([Rule 2](../xfoil-known-issues.md); `ci/validity-write-only.sh` enforces it).

**C cannot be folded in, and that is not an oversight.** Conditioning is a property of the problem
near that point, not of the arithmetic that was performed — there is nothing in a single run that
distinguishes a well-conditioned answer from one balanced on a pole. Seeing it requires perturbing
the input and running again: rerunning the reference with its panels jogged by one ULP (its
*twins*, `docs/conventions/terminology.md`). So a Class C point is reported `Valid`, because nothing about the
computation was invalid, and its uncertainty is carried by that study rather than by this flag.
Such a point is *ill-conditioned* in the sense of
[the project terminology](../conventions/terminology.md).

Where each is embodied:

| Class | Recorded in | Surfaces as a reason | Withholds numbers |
|---|---|---|---|
| A | `karman_tsien_margin_forces`, `karman_tsien_margin_pressure` | `KarmanTsienOutOfDomain` | yes |
| B | `mach_cl_newton_exhausted`, `converged`, `inviscid_cl_iterations` | `MachClNewtonExhausted`, `ClNewtonExhausted`, `ViscousNotConverged` | yes |
| C | — not in the record — | — | no |

On the [worked case](#worked-case) A and B between them separate the 61 points exactly, with no
false positives and no false negatives, and **both are load-bearing**: two points have a positive
denominator and are caught only by B. That is a measured result on one case, not a guarantee.

`SequenceHalted` belongs to none of the three. It is not a failure of a computation — it records
that no computation happened.

## Reason codes

A point is invalid on **exact, single-run facts only** — Class A or Class B above. Nothing is
inferred from the size or the smoothness of a returned number; see
[the caveat below](#what-valid-does-not-mean).

| `reason` | Class | What happened | What it means for the numbers |
|---|---|---|---|
| `KarmanTsienOutOfDomain` | A | The Kármán–Tsien denominator `β + BFAC·Cp_inc` reached zero or below at some node. The compressibility correction has a pole at `q/Q∞ = √(1 + 2β(1+β)/M∞²)`. | Past the pole Cp has **changed sign**, not merely lost accuracy. CL, CM and CDP are wrong, frequently in sign. Usually means the section is at too high an incidence for the Mach number, or that a CL-dependent type has pushed M∞ up to its 0.99 limit. |
| `MachClNewtonExhausted` | B | `SPECAL`'s CL(M) Newton used all 20 iterations without reaching `\|DCLM\| ≤ 1e-6`. | What is reported is iterate 20, not a solution. Only reachable under `TYPE 2`/`TYPE 3`, where the Mach number follows the lift. |
| `ClNewtonExhausted` | B | `SPECCL`'s alpha Newton used all 20 iterations without reaching `\|DALFA\| ≤ 1e-6`. | The alpha reported does not deliver the CL you asked for. Typically a CL beyond what the section can reach. |
| `ViscousNotConverged` | B | `VISCAL` finished without `LVCONV`: the viscous–inviscid iteration did not converge in `--max-iterations`. | The BL and the outer flow are not consistent with each other. Common past stall; raising `--max-iterations` sometimes helps and sometimes does not. |
| `SequenceHalted` | — | Only on `NotAttempted` points: the sweep stopped before this alpha. | Nothing was computed. |

## Diagnostics

`diagnostics` is present whenever the point was attempted, and is present **even when the numbers
are withheld** — these are facts about the computation, not physical results, and withholding them
would make the verdict unfalsifiable.

| Field | Meaning |
|---|---|
| `karman_tsien_margin_forces` | Smallest `β + BFAC·Cp_inc` over the force integration. **At or below zero is the violation**; a small positive value means the point is near the pole and hypersensitive. |
| `karman_tsien_margin_pressure` | The same over the stored Cp arrays. |
| `mach_cl_newton_exhausted` | `SPECAL`'s Newton hit its cap. |
| `cl_floored` | `MRCL` had a CL at or below zero and substituted `1e-6`: under a CL-dependent type there is no real Mach or Reynolds number for non-positive lift. |
| `mach_limited`, `re_limited` | `MRCL` limited the reported M∞ to 0.99, or Re to 100 × Re₁. |
| `converged`, `iterations`, `residual` | `LVCONV`, the VISCAL iteration count, and the final RMSBL. |
| `inviscid_cl_iterations` | `SPECCL`'s exit iteration (21 when exhausted). |

The margins are reported as numbers rather than collapsed to a flag on purpose. On the worked case
below, α = −1° reports CL = −52 and α = −5° reports CL = +0.128 — and the *first* is the smaller
violation of the two, spectacular only because its denominator is nearer zero. **The size of the
anomaly does not track the severity of the failure**, so any judgement based on how odd a number
looks will rank these backwards.

## The `cl_floored` field is not a validity test

`MRCL`'s flags describe the conditions the point was **reported at**, which is the last call
`SPECAL` makes — not every substitution made while iterating. Past the Kármán–Tsien pole the CL
comes back *positive* at negative incidence, so the final `MRCL` call is handed an ordinary-looking
number and clamps nothing. On the worked case, α = −30° records no substitution at all and is still
thoroughly out of domain. Read these fields as evidence, not as the verdict.

## What `Valid` does not mean

`Valid` means no formula was evaluated outside its domain (Class A) and no iteration ran out
(Class B). It does **not** mean the answer is well conditioned — **Class C is invisible to it**.

A converged, in-domain point can still sit where a 1-ULP change of input moves it by O(1) — next to
a pole, or at a station Newton that amplifies its seed. That cannot be detected from a single run at
all; it is measured with the twins, not by this flag. Post-CL~max~ polar points and
deep-stall solutions are the usual cases: often `Valid`, and often not reproducible to better than
their own noise. Treat validity as a *necessary* condition, not a sufficient one.

## Worked case

The NACA 64A010 swept −30° to +30° inviscid at `TYPE 2`, M₁ = 0.3, is the case the whole scheme was
built against, and it is written up in full at
[known issues §7.7](../xfoil-known-issues.md). `TYPE 2` makes the Mach number follow the lift, so
on a symmetric section the entire negative half of the sweep is outside the type's domain of
definition; `MRCL` substitutes a floored CL and a Mach limited to 0.99, where the pole sits at
`q/Q∞ = 1.153` — which this section's suction peak passes by α = ±1°.

Of the 61 points, 42 are invalid and 19 (α = +1° to +19°) are sound. XFOIL reports **positive lift
on a symmetric aerofoil at −30°** and prints nothing at all; its Newton fails outright at only six
of the 61. Both reason classes are load-bearing: α = −11° and α = 0° have a positive denominator and
are caught only by the exhaustion test.

## Checking yFoil's flagging against XFOIL

The instrumented reference build records the same occasions as branch events, so yFoil's flags can
be cross-checked against XFOIL directly — without comparing any numbers. The event names, and what
XFOIL does about each occasion on its own, are tabulated in
[the instrumentation guide](../xfoil-reference.md#validity-events).
