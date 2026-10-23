# Series cases

--8<-- "validation/series-cases/tables/run.md"

## What this study is for

These runs cover the aerofoil series and the flow regimes an aerodynamicist would normally use XFOIL
for, and show whether yFoil follows XFOIL through each of them. Every run is a full polar swept by
XFOIL's own polar procedure, from 0° up to +30°, then reinitialised and from 0° down to −30°. Each run
therefore goes past stall, which is where the two codes are hardest to keep together. Runs that exist
to show one of XFOIL's known issues are in [Known issues](../known-issues/index.md).

## How a run is compared

Each case is XFOIL's polar procedure at N = 160 panels and ITER 100: `ALFA 0 / ASEQ 1 30 1`, `INIT`,
`ALFA 0 / ASEQ -1 -30 -1`. That gives ITMAX + 5 iterations per point, and the sequence stops after
four consecutive unconverged points. An inviscid case is one SPECAL per alpha; a fixed-CL case is one
OPER `CL` command per point. yFoil is driven the same way, and every point XFOIL recorded is compared:

- **viscous points** with the studies' whole-run floor comparison (`studies/study-support/records.rs`).
  Iteration count, convergence, and the stagnation and transition stations must match exactly. Every
  per-iteration value and the converged point must be within `max(tol · scale, FLOOR_FACTOR · floor)`,
  where *floor* is how far XFOIL's own result moves under a 1-ULP change to its input;
- **inviscid points** on CL, CM, CDp and every node's γ, q and Cp, at `TOL_SOLVER`.

A point is *followed to the end* when every one of its iterations agrees with XFOIL's. *Parts at*
is the first point where yFoil's iterations leave XFOIL's. Whether that matters is decided on XFOIL
alone, using its **twins**: five reruns of XFOIL with every panel coordinate jogged by −1, 0 or +1 ULP
([terminology](../../development/conventions/terminology.md)).

- If the twins themselves spread by more than `DIVERGENCE_FLOOR` there, XFOIL is
  **ill-conditioned** at that point, and yFoil parting from it is expected. The result says nothing
  about the translation.
- If the twins agree, the parting is a **mismatch**. It is settled by replaying that one step: yFoil is
  seeded with XFOIL's exact state entering the call and must reproduce XFOIL's state after it
  (`tests/execution/steps.rs`).

!!! info "Studies show; tests check"

    This study demonstrates behaviour over whole sweeps. The tests check the same behaviour minimally,
    one step at a time ([testing conventions](../../development/conventions/testing.md)), and read
    nothing this study produces.

## Main results

!!! note "Stub — to be written"

    The headline across all runs: how far into each sweep yFoil follows XFOIL, and how many of the
    places where they part are XFOIL's own ill-conditioning.

--8<-- "validation/series-cases/tables/cases.md"

## Effect of accumulated error on run parting

A run [parts](../../development/conventions/terminology.md#run-parting) where yFoil stops matching
XFOIL, either because a branch decision differs or because a value leaves tolerance. Accumulated error
is not what causes it here: up to the point where a run parts, yFoil's final state stays within a small
fraction of its tolerance at every point, so the error does not build up from point to point, and the
parting happens within a single point, where XFOIL is ill-conditioned.

The table shows, for every run that keeps XFOIL's per-call state records, each point where yFoil parts
and the four points of the leg before it. *Parted by* says how the parting point parted (`branch: …` or
`value: …`). For each row the table gives the worst array of yFoil's final state
against XFOIL's: its worst entry (node number for a node array, `BL(is, ibl)` row for a station array,
side 1 first), the difference there, the twins' 1-ULP floor for that array, and the difference as a
ratio of its tolerance. A ratio well below 1 means the next point started from a state XFOIL's own
twins could not tell apart from XFOIL's.

--8<-- "validation/series-cases/tables/final-state.md"

## Reading the figures

There is one figure per run, with five rows, all against α:

1. CL and CD.
2. The last iteration's RMSBL, and the number of iterations taken.
3. Transition location x/c on the upper and lower surfaces. An open triangle marks a point where the
   final stagnation or transition station differs between the codes.
4. Where transition sits inside its station interval (0 at station ITRAN − 1, 1 at ITRAN).
5. yFoil's final state against XFOIL's, array by array: the worst array's difference as a ratio of its
   tolerance (left; the dotted line is the tolerance), and the difference against the 1-ULP floor
   (right). A state that stays at the floor from point to point has accumulated no error.

**XFOIL is red** (dashed, crosses) and **yFoil blue** (solid line, circles), so any visible
separation of the two colours is a difference between the codes. yFoil's circle is filled where every
iteration of the point agrees with XFOIL's. It is hollow where the iterations part within the point, or
where an earlier point in the same leg parted (that point started from a state XFOIL never computed,
so its comparison is informational).

- An **orange square** marks a point where that code went through a non-finite value. An orange plus
  at the top of a panel marks a value that is itself NaN and cannot be plotted.
- A **black diamond** marks where yFoil parts from XFOIL and XFOIL is ill-conditioned there (one of its
  own twins parts too).
- A **grey diamond** marks where yFoil parts from XFOIL and XFOIL is well-conditioned (a mismatch).
- Where the red line stops short of the blue one, XFOIL's sequence halted or hung.

## In this section

| Page | Runs |
|---|---|
| [Baselines](baselines.md) | NACA 0012 and NACA 4412 at Re 1e6 |
| [Reynolds number](reynolds-number.md) | Re 1e4 (with DAMP), 3e6, 1e7 and 1e8 |
| [Transition](transition.md) | Ncrit 4 and 11, and trips at x/c 0.22 |
| [Compressibility](compressibility.md) | M 0.5 and 0.6, including the NACA 16-212 |
| [Fixed CL and inviscid](fixed-cl-and-inviscid.md) | A fixed-CL sweep and an inviscid Kármán–Trefftz sweep |
