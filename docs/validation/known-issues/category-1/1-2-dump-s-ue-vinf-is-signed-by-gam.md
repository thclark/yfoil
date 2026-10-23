# 1.2 DUMP's `Ue/Vinf` is signed by GAM — Overcome (output only)

DUMP forms `UE = (GAM(I)/QINF)(1-TKLAM)/(1-TKLAM(GAM(I)/QINF)²)` (`xoper.f:1987`), so the lower surface
prints negative Ue. yFoil's `ue` column applies the same Kármán–Tsien transformation to UEDG (BLPRV's
`U2/QINF`), which is unsigned on both sides. `|DUMP Ue| == ue` on the surface.
