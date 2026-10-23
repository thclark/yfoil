# 1.1 DUMP prints one-step-stale closure quantities (H*, Cf, K, tau, Di) — Overcome (output only)

Verified in the source, not inferred:

- SETBL calls MRCHDU at the top of **every** Newton iteration (`xbl.f:93`, matched at
  `src/solver/setbl.rs:125`). MRCHDU re-marches and stores `THET, DSTR, UEDG, CTAU, MASS, TAU, DIS,
  CTQ, DELT, TSTR` (`xbl.f:1157-1167`). SETBL then evaluates the closures on that state
  (BLPRV/BLKIN/BLVAR, `xbl.f:219-220, 470`) and stores `TAU, DIS, CTQ, DELT, USLP` (`xbl.f:277-282`,
  matched at `setbl.rs:316-320`). Nothing writes TSTR in SETBL. Only after that do BLSOLV and UPDATE
  apply the Newton correction to `THET, DSTR, UEDG, CTAU, MASS`.
- So at the end of VISCAL the primaries are post-final-correction, while `TAU, DIS, CTQ, DELT, USLP,
  TSTR` describe the state *entering* the final iteration: they lag by exactly the last Newton
  correction.
- **Consequence for XFOIL's own output.** DUMP (`xoper.f:1955-1995`) prints `Ue, Dstar, Theta, H, HK,
  P, m` from current values but `Cf = TAU/(½Q²)`, `H* = TSTR/THET`, `K = TSTR·Ue³`, `tau`, `Di` from
  the lagged arrays. `H*` and `K` are the worst case: a lagged numerator over a current denominator.
  VPLO's `CF`, `CD` and `DELT` plots use `TAU, DIS, DELT` (`blplot.f:399, 499, 573`) and are lagged
  too; `H` is recomputed live with HKIN (`blplot.f:1385`), and the `TSTR/THET` version of that plot is
  commented out (`blplot.f:1387`), consistent with Drela knowing. `DT, DB, UE, N, CT` are live.
- Magnitude: the final correction, bounded by the convergence test (RMSBL < 1e-4 on the normalised
  corrections), so ~1e-4 relative or better on a converged point; unbounded on an unconverged one.
  On the reference case (NACA 0012, N=60, α=2°, Re=1e6) live and DUMP `H*` differ by ~1e-8.

**yFoil.** The stored arrays are reproduced exactly (`tests/subroutine/mrchdu.rs`,
`tests/execution/update.rs`, the SETBL tests); the solver is untouched. In the analysis JSON
(`BlSideOutput`, `src/output/foil.rs`) the canonical columns `hs, cf, cdis, delta, ctq, uslp` are
computed *live* by running XFOIL's own BLPRV → BLKIN → BLVAR on the converged primaries, and the lagged
arrays are emitted verbatim under `stored` (with `hs_dump`, `cf_dump` as DUMP prints them). A
side-by-side of live `hs` against DUMP's `H*` will differ at the ~RMSBL level by construction;
`tests/application/output.rs` bounds it. Column table: `docs/user-guide/reference/xfoil-reference.md`,
"Output" section.
