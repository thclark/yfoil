---
icon: lucide/shield-check
---

# Solution validity

XFOIL will hand back numbers for an operating point it has not actually solved. Sometimes it says
something on the console first; sometimes it says nothing at all. Either way the numbers arrive
looking exactly like good ones, and a polar plotted through them looks like a polar.

yFoil reproduces XFOIL's arithmetic exactly, including the cases above — that is the point of the
project ([Rule 2](known-xfoil-issues.md)) — so it computes the same wrong numbers. What it does
differently is **say so**, and withhold them unless you ask for them anyway. It also checks what
XFOIL cannot see from one run: whether the answer survives a 1-ULP jitter of the panel coordinates
(its [twins](#twins)).

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
  "classes": ["out-of-domain", "unconverged"],
  "diagnostics": { "karman_tsien_margin_forces": -0.0296, "mach_limited": true, ... }
}
```

`values` — the lift, drag and moment — is present only when the status is `Valid`, so it is absent
here. Pass
`--allow-invalid` to `yfoil analyse` or `yfoil polar` and it is filled in for invalid points too;
that flag changes **only** whether the numbers are present, never the status or the reasons, so a
result carrying numbers can still be recognised as one that should not be trusted.

For `yfoil analyse` the same rule covers `surface` (the q and Cp distributions) and
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
| `Valid` | Solved, every formula was evaluated inside its domain, every iteration met its tolerance, and the point's twins agree with it. `values` present. |
| `Invalid` | Solved, but at least one `reason` applies, so the numbers are not to be believed. `values` withheld unless `--allow-invalid`. |
| `NotAttempted` | Never solved. The sweep halted before reaching this alpha (by default XFOIL's NSEQEX rule, which ends a sequence after four consecutive non-converged points; see [halting](#halting-a-polar)). It is reported because you asked for it. |

## How failures are classified

Three distinct things can be wrong with a solved point, and they differ in something more important
than severity: **what it would take to know**. Keeping them apart is what stops the scheme from
either over-claiming or quietly ignoring a real problem.

| | Class | What it is | Provable from | Shape |
|---|---|---|---|---|
| **A** | Domain violation | A formula was evaluated outside its analytic domain, so its result is wrong rather than imprecise | a single run | binary, local, with a continuous margin |
| **B** | Iteration exhaustion | An iteration hit its cap without meeting its tolerance, so the value returned is not a solution to the problem posed | a single run | binary, global to the point |
| **C** | Conditioning | The point is converged and in domain, but sits where a 1-ULP jitter of input moves it by O(1) | **two or more runs**: the point and its twins | a spread, judged against two thresholds: a warning and an error |

**A and B** are exact facts about the computation as it was performed, so they are recorded as it
runs. They are detected by observation only: the solver records the evidence and never reads it
back, so no numeric result can depend on a validity classification
([Rule 2](known-xfoil-issues.md); `ci/validity-write-only.sh` enforces it).

**C cannot be seen in a single run.** Conditioning is a property of the problem near that point,
not of the arithmetic that was performed — there is nothing in one run that distinguishes a
well-conditioned answer from one balanced on a knife edge. Seeing it takes more runs with the input
perturbed, which is what the [twins](#twins) are: by default yFoil solves every point five more
times on jittered panels and compares. A point whose twins disagree is *ill-conditioned* in the
sense of [the project terminology](../../development/conventions/terminology.md), and its numbers
are withheld like those of an A or B point. Like A and B, the check observes and never steers: the
point reported is the point XFOIL computes, whatever its twins do.

C is judged against thresholds, where A and B are exact — but it is not a judgement on how odd a
number looks. The spreads fall into two populations separated by orders of magnitude: the noise
floor, where every twin reaches the point's own solution, and twins that reach another solution.
The error threshold sits in the gap between them, where moving it changes nothing; the warning
threshold sits in the upper reaches of the noise floor, and marks valid points whose twins agree
less closely than most ([the conditioning study](../../validation/conditioning/index.md)).

Where each is embodied:

| Class | Recorded in | Surfaces as a reason | Withholds numbers |
|---|---|---|---|
| A | `karman_tsien_margin_forces`, `karman_tsien_margin_pressure` | `KarmanTsienOutOfDomain` | yes |
| B | `mach_cl_newton_exhausted`, `converged`, `inviscid_cl_iterations`, `nonfinite_station_failures` | `MachClNewtonExhausted`, `ClNewtonExhausted`, `ViscousNotConverged`, `NonFiniteStationFailure` | yes |
| C | `conditioning` (the twins' spreads) | `IllConditioned` (beyond the error threshold; beyond only the warning threshold, an `IllConditioned` warning) | yes, unless `--conditioning off` |

On the [worked case](#worked-case) A and B between them separate the 61 points exactly, with no
false positives and no false negatives, and **both are load-bearing**: two points have a positive
denominator and are caught only by B. That is a measured result on one case, not a guarantee.

`SequenceHalted` belongs to none of the three. It is not a failure of a computation — it records
that no computation happened.

A point already invalid for an A or B reason is not judged on its twins: there is nothing left to
protect, and its record carries no `conditioning`.

## Reason codes

A point is invalid on exact, single-run facts — Class A or Class B above — or because its twins
disagree with it (Class C). Nothing is inferred from the size or the smoothness of the point's own
numbers; see [the diagnostics](#diagnostics).

| `reason` | Class | What happened | What it means for the numbers |
|---|---|---|---|
| `KarmanTsienOutOfDomain` | A | The Kármán–Tsien denominator `β + BFAC·Cp_inc` reached zero or below at some node. The compressibility correction has a pole at `q/Q∞ = √(1 + 2β(1+β)/M∞²)`. | Past the pole Cp has **changed sign**, not merely lost accuracy. CL, CM and CDP are wrong, frequently in sign. Usually means the section is at too high an incidence for the Mach number, or that a CL-dependent type has pushed M∞ up to its 0.99 limit. |
| `MachClNewtonExhausted` | B | `SPECAL`'s CL(M) Newton used all 20 iterations without reaching `\|DCLM\| ≤ 1e-6`. | What is reported is iterate 20, not a solution. Only reachable under `TYPE 2`/`TYPE 3`, where the Mach number follows the lift. |
| `ClNewtonExhausted` | B | `SPECCL`'s alpha Newton used all 20 iterations without reaching `\|DALFA\| ≤ 1e-6`. | The alpha reported does not deliver the CL you asked for. Typically a CL beyond what the section can reach. |
| `ViscousNotConverged` | B | `VISCAL` finished without `LVCONV`: the viscous–inviscid iteration did not converge in `--max-iterations`. | The BL and the outer flow are not consistent with each other. Common past stall; raising `--max-iterations` sometimes helps and sometimes does not. |
| `NonFiniteStationFailure` | B | A station Newton failed with a non-finite residual, and the march was carried over it by XFOIL's iterate fallback: beyond the third station the station is set from the one upstream, at the first three the non-finite iterate is kept. | The run passed through a NaN and kept going. See [known issues §4](known-xfoil-issues.md), the TRCHEK2 row. |
| `IllConditioned` | C | A 1-ULP jitter of the panel coordinates moves a reported quantity (CL, CD, CM, …) further than the [error threshold](#twins): at least one twin reached another solution. | The point is balanced on a knife edge: another solution is one rounding error away, so the one reported is an accident of the last bit of the input. Usual past stall and in deep separation. `--conditioning off` reports it as XFOIL would. |
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
| `nonfinite_station_failures` | Station Newton failures with a non-finite residual. |
| `conditioning` | How far the [twins](#twins) moved what the point reports: `seeds`, `thresholds` (`warning` and `error`), `spread` (one entry per reported quantity), `worst` (the quantity with the largest spread) and `twins_not_valid`. Absent with `--conditioning off`, or when the point is invalid for an A or B reason. |

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

## Twins

A **twin** of a point is the same solve — the same command, from the same initial state — on panels
whose every coordinate has been jittered by −1, 0 or +1 ULP, x and y independently, from a fixed
seed. yFoil solves each point with five twins (seeds 1 to 5, the same jitter the validation studies
apply to XFOIL) unless `--conditioning off` is given, and compares.

- **The spread** of a reported quantity is the furthest any twin moved it, `|a − b| / max(|a|, |b|, 1)`:
  absolute for the coefficients, relative for anything larger. A twin that went non-finite where the
  point did not has spread `null`, the largest there is.
- **The point is ill-conditioned** when any spread exceeds the error threshold, 1e-4 by default
  (`--spread-error-threshold`, `CONDITIONING_ERROR_THRESHOLD`). It is then `Invalid` with
  `IllConditioned`, and its numbers are withheld. Twins that agree mostly do so to the noise floor,
  around 1e-12; a twin that disagrees has reached another solution, and differs by 1e-2 or more.
- **The point may be ill-conditioned** when its largest spread exceeds only the warning threshold,
  1e-8 by default (`--spread-warning-threshold`, `CONDITIONING_WARNING_THRESHOLD`). It stays `Valid`
  with its numbers, and carries a warning naming the quantity and its spread —
  `"warnings": [{"IllConditioned": {"quantity": "cl", "spread": 1.86e-7}}]`, printed as *solution
  may be ill-conditioned: its twins spread cl by 1.86e-7*. Such points reached the same solution,
  but less closely than most: typically a twin's convergence test tipped by an iteration.
- **Where the thresholds sit**, measured over the series and branch-coverage cases, is
  [the conditioning study](../../validation/conditioning/index.md). Both can be set for one run; the
  warning threshold must be no larger than the error threshold.
- **Only values decide.** A twin that needed more iterations, or did not converge at all, is
  compared by the values it ended with. When the point converged and some twins did not, the record
  carries a warning — `"warnings": [{"TwinsNotValid": {"count": 2}}]`, printed as *solution converged
  but 2 twins did not* — which changes nothing about the status but lets you filter on it.
- **Per point, from the point's own initial state.** In a polar a point starts from its
  predecessor's converged boundary layer, and so do its twins: they are twins of that one solve, not
  of the whole sweep. A point's verdict therefore depends only on what its `id` hashes — its inputs
  and its predecessor — so it is the same however the point was reached, and results from different
  runs can be merged.
- **Twins observe; they never steer.** The numbers reported are the base solve's, identical with
  or without twins. By default the sweep's halting rule counts non-converged points only; a polar
  stops on ill-conditioned points only when asked to ([halting](#halting-a-polar)).

Twins cost five more solves per point, run in parallel: about twice the wall-clock time on a
machine with six or more cores. `--conditioning` sets what they do:

| `--conditioning` | Twins run | Ill-conditioned points | Twins' records |
|---|---|---|---|
| `on` (default) | yes | `Invalid (IllConditioned)` | — |
| `off` | no; no `conditioning` is recorded, and every converged, in-domain point is reported as XFOIL would report it | — | — |
| `study` | yes | `Invalid (IllConditioned)` | in a `study` array |

A **study** keeps everything the check computed. Each twin is written as a record of its own beside
the point's: `role: twin`, `twin_of` (the point's `id`), its `seed`, and an `id` of its own (a twin
hashes its jittered panels and seed, and cites the point's predecessor, the state it started from).
In a polar the `study` array lists every point's twins; for `yfoil analyse -o` and the records
`yfoil polar --distributions` embeds, each record's `study` holds its twins' full records, geometry,
surface and boundary layer included. **A study's records always carry their numbers**, whatever
their status and whatever `--allow-invalid` says: withholding them would hide exactly what a study
is for, and each still carries its own `status` and `reasons`. `results` is unchanged by a study. [Polars](../commands/polars.md#keeping-the-twins-studies) works through an
example, and [halting a leg](../commands/polars.md#halting-a-leg) shows a polar stopped at its first
ill-conditioned point.

`--allow-invalid` is separate from all three: with it the twins still run and an ill-conditioned
point is still `Invalid`, but its numbers are shown in `results`.

## Classes on records

Every record lists the classes of its reasons in `classes`, once each in the order A, B, C, by the
names the command line uses:

| Class | Name | Reasons |
|---|---|---|
| A | `out-of-domain` | `KarmanTsienOutOfDomain` |
| B | `unconverged` | `MachClNewtonExhausted`, `ClNewtonExhausted`, `ViscousNotConverged`, `NonFiniteStationFailure` |
| C | `ill-conditioned` | `IllConditioned` |

## Halting a polar

`yfoil polar --halt-on CRITERION[=N],...` halts a leg of the sweep after N consecutive points
that meet a criterion: a class name (N defaults to 1), `xfoil` (N defaults to 4), or `none`. Each
criterion keeps its own count, reset by a point that does not meet it. The default, `xfoil`, is
XFOIL's own NSEQEX rule — four consecutive viscous points whose iteration did not converge — and is
narrower than `unconverged`, which also counts the inviscid Newtons and a converged point that
passed through a NaN; only `xfoil` reproduces XFOIL's sweep exactly. Halting decides which alphas
are attempted, never the numbers of a point that is: the alphas past a halt are `NotAttempted`, and
the polar's `halts` records each halted leg's criterion and its last alpha.

A fourth class, D (`anomalous`: a point inconsistent with its neighbours, known only from
neighbouring operating points), is designed in
[issue #12](https://github.com/thclark/yfoil/issues/12) and not yet available.

## What `Valid` does not mean

`Valid` means no formula was evaluated outside its domain (Class A), no iteration ran out
(Class B), and none of five twins moved the answer beyond the error threshold (Class C). A valid
point with an `IllConditioned` warning passed that test with less room than most. Five twins are a
sample, not a proof: a point can be closer to a knife edge than they happened to reach, and a point
solved with `--conditioning off` has not been checked at all. Treat validity as strong evidence, not a
guarantee.

## Worked case

The NACA 64A010 swept −30° to +30° inviscid at `TYPE 2`, M₁ = 0.3, is the case the whole scheme was
built against, and it is written up in full at
[known issues §7.7](known-xfoil-issues.md). `TYPE 2` makes the Mach number follow the lift, so
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
[the instrumentation guide](xfoil-reference.md#validity-events).
