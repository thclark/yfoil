# 6.5 TGAP on a sharp trailing edge does not deliver the requested gap — Replicated

`xgdes.f:1238-1244` (TGAP): when the buffer airfoil's trailing-edge gap is zero, the direction
along which the two surfaces are moved apart is `(DXU, DYU) = (−½(YBP(NB) − YBP(1)), ½(XBP(NB) −
XBP(1)))`, the mean of the two end tangents of the spline. That vector is unit length only for a
cusp: for a trailing-edge angle τ it has length cos(τ/2) (times the small stretch of the
spline-parameter tangents), and the surfaces end up `cos(τ/2)` × the requested gap apart. With a
gap already open the direction is the existing gap's unit vector and the requested gap is
delivered exactly. Measured on the tracked `naca63-415_n160_tgap` case (τ ≈ 7°): `TGAP 0.002`
gives 2.0012e-3. yFoil's `set_te_gap` reproduces the moved nodes bitwise
(`tests/subroutine/tgap.rs`); the resulting gap is asserted against the requested one only when a
gap existed. The gap TGAP produces is what `TECALC` then measures, so the `SHARP` decision and
the base-drag treatment see the delivered gap, not the requested one.
