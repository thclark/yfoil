# Physical invariants

Two codes can share a misunderstanding: a translation that reproduces XFOIL exactly reproduces any
error XFOIL makes. So yFoil is also checked against properties that hold whatever XFOIL does — of
the geometry, the panel method and the boundary-layer closures. These are the `invariants` tests
(`tests/invariants/physics.rs`); they need no fixture and compare with no reference.

## What is checked

| Property | What the test requires |
| --- | --- |
| A symmetric section is symmetric | the NACA 0015's upper and lower extremes are equal and opposite, to 1e-3 |
| A cambered section has positive camber | the NACA 6412 and 23018 reach further above the chord line than below it |
| A symmetric section carries no lift at zero incidence | inviscid \(C_L\) of the NACA 0012 at α = 0 is zero to 1e-12 |
| A cambered section lifts at zero incidence | inviscid \(C_L\) of the NACA 4412 at α = 0 exceeds 0.2 (thin-aerofoil theory gives about 0.5) |
| The lift-curve slope is near thin-aerofoil theory's 2π | the inviscid NACA 0012's \(\mathrm{d}C_L/\mathrm{d}\alpha\) between 0° and 5° lies between 5 and 7 per radian, sharp and blunt trailing edge alike, the two within 15 % of each other |
| The Kutta condition holds | the trailing-edge vorticities of the unit solutions cancel to 1e-9, sharp and blunt trailing edge |
| A sharp trailing edge's vorticity is smooth | the curvatures of the vorticity at the two trailing-edge nodes are antisymmetric to 1e-9 |
| A blunt trailing edge gives sensible results | the inviscid NACA 0012 with a 0.2 % chord trailing-edge gap at α = 0 has \(|C_L| < 0.05\), \(|C_{D,p}| < 0.01\) and no surface speed above 3 |
| The laminar closures hold at the Blasius profile | at \(H = 2.591\): \(H_k = H\) at zero Mach to 1e-10, \(H^*\) between 1.5 and 1.7, \(C_f\,Re_\theta\) of order one at \(Re_\theta\) from 500 to 5000 |
| Amplification switches on at the critical Reynolds number | the amplification rate is zero below the critical \(Re_\theta\) and positive above it |

## What is not checked yet

CLAUDE.md's Rule 6 names three invariants that the tests above only approach:

- **CL = CM = 0 at α = 0 on a symmetric section, viscous, to the noise floor.** The inviscid lift
  is held to 1e-12; no viscous point is checked, and \(C_M\) not at all.
- **A mirrored section at −α gives the mirrored solution.** No test yet.
- **The Blasius flat plate.** The closures are checked at the Blasius shape factor; a marched
  flat-plate boundary layer is not.

Until they land, the strongest evidence that yFoil's physics is right remains its agreement with
XFOIL, the code whose physics is trusted:
[agreement with XFOIL under normal conditions](agreement-with-xfoil-under-normal-conditions.md).
