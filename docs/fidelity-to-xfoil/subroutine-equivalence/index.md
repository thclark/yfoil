# Subroutine equivalence

The smallest unit of the comparison: one translated XFOIL subroutine, run on the inputs XFOIL
handed it in a real run, its outputs required to match XFOIL's. Nothing upstream of the call is
yFoil's, so a difference is in that subroutine and nowhere else. These are the `subroutine` tests
(`tests/subroutine/`), each module naming the fixture it reads and the command that regenerates
it; the tolerances are the [named constants](../tolerances/index.md#named-tolerance-constants).

| Page | Routines | Test modules |
| --- | --- | --- |
| [Boundary-layer closures](boundary-layer-closures.md) | HKIN, CFL, HSL, DIL, HST, CFT (and AXSET, BLKIN, BLVAR in the tests) | `closures.rs` |
| [Transition and amplification](transition-and-amplification.md) | DAMPL, TRCHEK2 | `closures.rs`, `transition.rs` |
| [Boundary-layer marching and the Newton system](boundary-layer-marching-and-the-newton-system.md) | MRCHUE, MRCHDU, SETBL, BLSOLV, SPECCL | `mrchue.rs`, `mrchdu.rs`, `setbl.rs`, `blsolv.rs`, `speccl.rs` |
| [Inviscid panel routines](inviscid-panel-routines.md) | GGCALC, QDCALC, XYWAKE, SETEXP, QWCALC, the station pointers and the inviscid velocity layer | `ggcalc.rs`, `xywake.rs`, `pointers.rs` |
| [XFOIL's own geometry routines](xfoils-own-geometry-routines.md) | NACA4, NACA5, PANGEN, TGAP | `pangen.rs`, `tgap.rs`, `pane_legacy.rs` |

The first two pages are generated from the closure fixtures, with every call's error; the others
describe their tests, which assert rather than tabulate. A subroutine whose result depends on the
iterations before it — UPDATE, the VISCAL loop — is tested as a step of a run instead:
[execution equivalence](../execution-equivalence/index.md).
