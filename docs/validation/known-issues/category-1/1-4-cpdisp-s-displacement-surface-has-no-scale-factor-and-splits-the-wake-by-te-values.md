# 1.4 CPDISP's displacement surface has no scale factor and splits the wake δ* by TE values — Replicated, extended

`xplots.f:699-754`: the displacement surface is `X + NX·DSTR` at true geometric scale; the wake δ* is
split into upper/lower fractions `DSF1 = (DSTR(IBLTE(1),1) + ½ANTE)/DSTR(IBLTE(2)+1,2)` (and DSF2
likewise) with a 0.5/0.5 fallback when the first wake δ* is exactly zero, and the wake normal points
to the *lower* side so the upper edge is `X − N·DSTR·DSF1`. `src/output/foil_plot.rs` follows this
exactly for δ*, and additionally (yFoil's own): a per-quantity scale factor, the same construction for
other quantities with a 0.5/0.5 split, and markers for stagnation, transition and separation, which
XFOIL never draws on the airfoil plot. The separation/reattachment markers are derived from the sign of
the live Cf; XFOIL reports no separation location anywhere.
