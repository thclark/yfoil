| file | lines | misaligned | max abs diff | max rel (abs>1e-8) | p99 rel | median rel | stage |
|---|---|---|---|---|---|---|---|
| xfoil_panels.dat | 71 | 0 | 7.30e-15 | 5.25e-14 | 2.45e-14 | 1.45e-16 | geometry after ABCOPY |
| xfoil_inviscid.dat | 129 | 0 | 1.95e-12 | 4.12e-11 | 9.49e-12 | 5.06e-14 | inviscid GAM/QINV |
| xfoil_dij.dat | 5628 | 0 | 1.01e-08 | 5.89e-06 | 2.59e-07 | 1.64e-11 | DIJ matrix |
| blsolv_input.dat | 17802 | 0 | 6.49e-05 | 4.54e-03 | 1.36e-06 | 4.91e-11 | SETBL VA/VB/VM/VDEL (calls 1..3) |
| blsolv_output.dat | 231 | 0 | 1.25e-10 | 2.13e-08 | 1.55e-08 | 5.33e-11 | BLSOLV VDEL out (calls 1..3) |
| viscal_iter.dat | 120 | 0 | 1.27e-11 | 1.05e-10 | 5.38e-11 | 3.91e-13 | per-iteration CL/CD/RMSBL/RLX (VISCAL call 1) |
| viscal_final.dat | 169 | 0 | 9.47e-12 | 1.55e-10 | 1.14e-11 | 7.97e-14 | converged point (VISCAL call 1) |
| viscal_iters_all.dat | 6 | 0 | 2.50e-12 | 1.05e-10 | 3.19e-11 | 2.94e-13 | per-iteration record, all VISCAL calls |
| viscal_points.dat | 21 | 0 | 2.18e-13 | 5.41e-12 | 2.67e-12 | 1.00e-13 | per-point CL/CD/CM/XTR (all VISCAL calls) |
