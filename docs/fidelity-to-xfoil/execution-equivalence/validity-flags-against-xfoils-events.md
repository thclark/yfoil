# Validity flags against XFOIL's events

yFoil marks a point `Invalid` on exact, single-run facts: a formula evaluated outside its domain, or
an iteration that ran out ([failure classes A, B and C](../../using-yfoil/whether-to-trust-a-result.md#failure-classes-a-b-and-c)).
XFOIL computes the same points and mostly says nothing. This page is how yFoil's flags are checked
against the reference anyway — by the branches XFOIL took, not by its numbers.

## The branch events

`xfoil/instrumentation/instrument/19-validity-events.patch` records, as `EVLOG` branch events in
`events.dat`, every occasion on which XFOIL goes on to deliver numbers for a point it has not
properly solved. It exists so that yFoil's flagging can be checked against the reference **without
comparing any numbers** (the flags are [whether to trust a result](../../using-yfoil/whether-to-trust-a-result.md)).

The reason the events are needed at all is that XFOIL mostly does not record these occasions
anywhere a fixture can see. What it does about each, on its own:

| Event | What XFOIL does about it today |
|---|---|
| `CLCALC_KT_DOMAIN` | **nothing at all** — forms the denominator inline, never tests it |
| `MRCL_CL_FLOOR` | **nothing at all** — silently substitutes `CLA = MAX(CLS, 1e-6)` |
| `CPCALC_KT_DOMAIN` | console only (`Local speed too large`) |
| `SPECAL_MINF_FAIL` | console only (`Minf convergence failed`) |
| `MRCL_MACH_LIMIT` / `MRCL_RE_LIMIT` | console only |
| `SPECAL_ENTER` | — (marker, carries alpha) |

The console messages are not a substitute: `stdout.txt` is not tracked, and `specal_points.dat`
carries no `ITCL`, so for an inviscid case the event is the only record in a fixture that the CL(M)
Newton was exhausted.

`SPECAL_ENTER` matters more than it looks. Events carry `NSETBLC`, the SETBL call counter, which is
**zero throughout an inviscid run** — so without a marker an inviscid event could not be attributed
to an operating point. Every event between one marker and the next belongs to that SPECAL call, and
the marker's real value is the alpha.

## The check, measured

Measured on `naca64a010_n240_inviscid_polar30_type2_m03` (2026-09-17), event counts against the
console messages they correspond to:

```
 360 CLCALC_KT_DOMAIN      <- no console equivalent
  80 CPCALC_KT_DOMAIN      = "Local speed too large"      80
 181 MRCL_CL_FLOOR         <- no console equivalent
 360 MRCL_MACH_LIMIT       = "limiting Mach"             360
  61 SPECAL_ENTER          = 61 alphas
   6 SPECAL_MINF_FAIL      = "Minf convergence failed"      6
```

Exact agreement everywhere XFOIL speaks. Comparing per alpha against yFoil's own flags on the same
case: the Newton-exhaustion sets are identical (`[-11, -8, -3, -2, -1, 0]`) and the verdicts are
identical at all 61 points (42 invalid, α = +1°…+19° sound).

## Event and flag differ by design

**Granularity differs from yFoil's record by design.** An event fires on *any* call within the span,
while yFoil's fields describe the *reported* state — the last call. So a point whose intermediate
Newton iterate left the domain but whose final one did not carries the event and not the flag. On
the case above that happens at exactly one point, α = −11°, which both codes flag invalid anyway
through the exhaustion test. A comparison written against these events should be at verdict level,
or should account for the difference.

As with the rest of the series, the patch is inert
([proof that instrumentation changes nothing](../the-reference-xfoil-build/proofs-that-the-reference-is-sound.md#proof-that-instrumentation-changes-nothing)).

## The test

`tests/known_issues/kt_pole_7_7.rs` checks yFoil's record against these events on
`naca4412_n60_inviscid_m07_a10` and `naca64a010_n60_inviscid_type2_m03_am23_a0_a19`, at verdict
level. The case the scheme was built against is
[known issues §7.7](../xfoils-known-issues/the-issue-register.md).
