# 6.3 CHORD and APANEL are not what their names suggest — Replicated

`CHORD` (GEOPAR) is the distance from the LE point *on the spline* to the TE midpoint, 0.99995721 for
the 60-panel NACA 0012 rather than 1; `APANEL` (APCALC) is the panel-*normal* angle, 3π/2 from a
node-tangent angle. Both were yFoil misreadings found by the S3 tests; both are now XFOIL's, and
`CHORD` is what SETEXP scales the wake by.
