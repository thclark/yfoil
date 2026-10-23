# 2.3 TRDIF's transition-point station reuses station-2 COMMON — Replicated

`xblsys.f:355` "temporarily set '2' variables from 'T' for BLKIN": XFOIL overwrites
`X2/T2/D2/U2/AMPL2/S2` in place without calling BLPRV, so `U2_UEI`, `U2_MS` and `DW2` at the transition
point are still station 2's (`src/bl/system.rs:2191-2193`).
