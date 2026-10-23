# 6.6 NACA generator constants: the ordinate program's precision — Reference note

Not XFOIL: the NASA/PDAS `naca456` program that the geometry generators are tested against
(`tests/fixtures/naca456/`) carries `PI = 3.141592654` in its 6-series mean lines and
`3.14159265` in the φ grid of the 6-series mapping, and inverts its arc-length spline with Brent's
method at a 1e-6 tolerance, reporting the ordinate at the station reached. yFoil uses π and inverts
to round-off; the tolerances in `tests/common/utilities/tolerances.rs` are derived from those three facts.
naca456 also reports the aft slope of the 4-digit modified thickness form as dy/d(1 − x).
