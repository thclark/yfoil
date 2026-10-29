# Conditioning

--8<-- "validation/conditioning/tables/run.md"

!!! info "Used by yFoil"

    The two thresholds every analysed point is judged against — `CONDITIONING_WARNING_THRESHOLD`
    (1e-8) and `CONDITIONING_ERROR_THRESHOLD` (1e-4), in `src/output/conditioning.rs` — were set by
    hand from this study's findings. Nothing reads the study's output, so rerunning it changes
    nothing until the constants are revisited by hand. Both can be set for one run
    (`--spread-warning-threshold`, `--spread-error-threshold`); what they mean for a result is in
    [Solution validity](../../user-guide/reference/solution-validity.md#twins).

## What this study is for

A converged, in-domain point can still be balanced on a knife edge: another solution lies one
rounding error away, and the one reported is an accident of the last bit of the input. One run
cannot show this. yFoil therefore solves each point five more times by default, on copies of the
panels given a 1-ULP jitter (its **twins**), and compares. Comparing needs a rule: how far may a
twin move a reported quantity before the point is no longer to be believed?

This study answers that question, and only that one. It measures how far the twins move every
point of a broad set of operating points, and shows where thresholds can sit so that the verdict
depends on the points and not on the exact threshold chosen. It uses yFoil alone; nothing the
reference computed is read.

## How

**The points.** Two groups of cases from `xtask/fixtures-config/cases.toml` supply panels and an
OPER procedure, which yFoil runs with its twins:

- the **series** cases (`group = "series"`), the aerofoil series and flow regimes of interest under
  normal conditions, each swept 0 → ±30° by XFOIL's polar procedure (see
  [Series cases](../series-cases/index.md));
- the **branch-coverage** cases (`group = "branch-coverage"`), single points and short sequences
  that take XFOIL's rarely-taken branches: inviscid, fixed CL, `TYPE 2` and `3`, forced transition,
  `DAMP` (see [Branch coverage](../branch-coverage/index.md)).

**The measure.** Let $q_v$ be a reported quantity $v$ of the point and $q_{v,k}$ the same
quantity from its twin $k = 1, \dots, 5$. The *spread* of $v$ is the furthest any twin moved it,
and the point's *largest spread* is the largest over its quantities:

$$
s_v = \max_k \frac{|q_v - q_{v,k}|}{\max\left(|q_v|,\ |q_{v,k}|,\ 1\right)},
\qquad
S = \max_v s_v .
$$

The quantities $v$ are α, CL, CM, CDp, CD, CDf and the x/c of the upper and lower transition
(the first four only for an inviscid point). Each twin is compared with the point, never with
another twin. The floor of 1 in the denominator makes the spread an absolute difference for every
quantity of order one or less, which is all of them in practice, and a relative one only beyond
that. A spread that is not a number (a twin went through a NaN where the point did not) counts as
the largest possible, so then $S$ is too. This is the measure yFoil applies to every point
(`src/output/conditioning.rs`); $S$ is what the figures below plot.

**Which points are measured.** A point already invalid on its own run, because it did not converge
or left the Kármán–Tsien domain, is withheld whatever its twins do. Its twins are not compared, so
it is counted but not measured.

--8<-- "validation/conditioning/tables/cases.md"

## Main results

**The spreads fall into two populations.** Most points' twins agree to the noise floor of a
converged Newton solve: at or below 1e-10, the floor the [tolerances study](../tolerances/index.md)
measures for converged state. A small number agree less closely, by up to 1e-5. These reached the same solution, but some
twin's convergence test tipped by an iteration on the way. Then comes an empty band of decades. Above it are the
points where a twin reached *another* solution, moving a coefficient by 1e-2 or more, or making it
non-finite.

<ironlab-figure src="figures/spread-histogram.fig" format="fig" name="spread-histogram" alt="Points per decade of their largest spread over five 1-ULP twins, series and branch-coverage cases, with the warning and error thresholds">Points per decade of their largest spread over five 1-ULP twins, series and branch-coverage cases, with the warning and error thresholds</ironlab-figure>

*The number of measured points in each decade of their largest spread: series cases as blue
circles, branch-coverage cases as green squares. A spread that is not a number is drawn in orange
at 1e2. The dash-dotted line is the warning threshold, the dashed line the error threshold.*

--8<-- "validation/conditioning/tables/histogram.md"

**The verdict does not depend on where in a gap a threshold sits.** Counting the points beyond each
threshold from 1e-14 to 1 gives a curve with flat steps. On a step, the same points are beyond the
threshold however it moves. Each of yFoil's thresholds sits on one.

<ironlab-figure src="figures/threshold-plateau.fig" format="fig" name="threshold-plateau" alt="Points whose largest spread exceeds each threshold from 1e-14 to 1, with the warning and error thresholds and the gaps around them">Points whose largest spread exceeds each threshold from 1e-14 to 1, with the warning and error thresholds and the gaps around them</ironlab-figure>

*The number of measured points beyond each threshold. The dotted lines are the edges of the gap
around each of yFoil's thresholds: the largest spread at or below it and the smallest above it.
Across each gap, the same points are beyond the threshold.*

--8<-- "validation/conditioning/tables/gap.md"

## The thresholds

**The error threshold, 1e-4**, decides whether a point is *ill-conditioned*: `Invalid` with reason
`IllConditioned`, its numbers withheld. It sits in the gap between the two populations. Below it
the twins all reached the point's own solution; above it at least one reached another. It is set
nearer the noise floor than the gap's centre for two reasons. A twin that reaches a neighbouring
solution need not move a coefficient far. And a noise-floor point would need a spread
several times larger than any measured here to cross it.

**The warning threshold, 1e-8**, marks valid points that *may be* ill-conditioned: they keep their
status and numbers, and carry an `IllConditioned` warning naming the quantity and its spread. It
separates the bulk of the noise floor from the points whose twins agree less closely than a
converged solve's floor. Those points reached the same solution, but some twin's convergence tipped
on the way. They are worth knowing about, not worth withholding, so a user can filter on the
warning.

**No quantity sets the error threshold alone.** For every reported quantity, the largest spread
among the points within the error threshold is below it. Every ill-conditioned point is flagged by
CL, CM, CDp or CD, so no verdict rests on the transition locations alone.

--8<-- "validation/conditioning/tables/gap-by-quantity.md"

## Where the flagged points are

Ill-conditioned points lie past stall, on the negative leg of the cambered sections' polars and
at the onset of separation in compressible flow. Points that may be ill-conditioned lie mostly at
transition jumps and deep in the post-stall leg. No branch-coverage point is ill-conditioned.

<ironlab-figure src="figures/spread-against-alpha.fig" format="fig" name="spread-against-alpha" alt="Every measured point's largest spread over five 1-ULP twins against alpha, series and branch-coverage cases, with the warning and error thresholds">Every measured point's largest spread over five 1-ULP twins against alpha, series and branch-coverage cases, with the warning and error thresholds</ironlab-figure>

*Every measured point's largest spread against α: series cases as blue circles, branch-coverage
cases as green squares. Open points are within the warning threshold, half-filled points between
the thresholds, filled points beyond the error threshold. Points invalid on their own run are grey
diamonds along the bottom.*

### Ill-conditioned (withheld)

--8<-- "validation/conditioning/tables/ill-conditioned.md"

### May be ill-conditioned (warning)

--8<-- "validation/conditioning/tables/may-be-ill-conditioned.md"

## The plateau in full

--8<-- "validation/conditioning/tables/plateau.md"
