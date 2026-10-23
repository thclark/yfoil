# 7.6 Single precision and FPE checks — Reference-build facts

Stock XFOIL is single precision and agrees with the double-precision reference only to ~1e-7. The
gfortran Makefile ships with the bounds/FPE checks commented out (`bin/Makefile_gfortran:120`;
`bin/Makefile` has them on). The DP reference built with `-finit-real=snan -ffpe-trap=invalid` runs the
smoke case clean; the only raised (untrapped) flag is `IEEE_DIVIDE_BY_ZERO` from `PLTINI`
(`xplots.f:37`), plot-scale setup that runs even with graphics off, outside the solver.
