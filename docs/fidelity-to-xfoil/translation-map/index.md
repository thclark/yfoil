---
icon: lucide/library
---

# Translation map

Where each XFOIL variable and routine lives in yFoil, and how XFOIL itself is put together. This
section is for reading the two codes side by side: when a fixture disagrees, it is how you find
the Rust that translates the Fortran in question, and the XFOIL state the Rust mirrors.

## XFOIL in one paragraph

XFOIL (created by Mark Drela at MIT, 1980s) performs viscous-inviscid analysis of 2D airfoil sections using:

- Panel method with vortex distribution for inviscid flow
- Integral boundary layer solver with e^N transition prediction
- Coupled viscous-inviscid iteration (VISCAL)

## XFOIL's source files

| File       | Purpose                      | Key Subroutines                     |
|------------|------------------------------|-------------------------------------|
| `xoper.f`  | Main operations and coupling | VISCAL, SPECAL, SPECCL              |
| `xblsys.f` | BL Newton system             | BLPRV, BLKIN, BLVAR, BLSYS, TRCHEK2 |
| `xbl.f`    | BL marching                  | SETBL, MRCHUE, MRCHDU, UPDATE       |
| `xpanel.f` | Panel method                 | PSILIN, QDCALC, UICALC, QVFUE       |
| `spline.f` | Cubic splines                | SPLINE, SEVAL, SPLIND               |
| `xfoil.f`  | Panelling                    | PANGEN, GETPAN (PPAR), TECALC       |
| `xgeom.f`  | Geometry                     | LEFIND, SCALC, SEGSPL               |
| `naca.f`   | NACA generation              | NACA4, NACA5                        |

## The main routines in yFoil

The full XFOIL → yFoil name mapping is [variable and routine names](variable-and-routine-names.md);
the naming rules are [naming](../../contributing/conventions/naming.md). Where the main
translated routines live:

| XFOIL             | yFoil                                                  |
|-------------------|--------------------------------------------------------|
| `xoper.f:VISCAL`  | `src/solver/viscal.rs::solve_viscous`                  |
| `xoper.f:SPECAL`  | `src/solver/specal.rs::solve_inviscid_at_alpha`        |
| `xbl.f:SETBL`     | `src/solver/setbl.rs::assemble_newton_system`          |
| `xbl.f:UPDATE`    | `src/solver/update.rs::apply_newton_update`            |
| `xbl.f:MRCHUE`    | `src/bl/mrchue.rs::march_direct`                       |
| `xbl.f:MRCHDU`    | `src/bl/mrchdu.rs::march_prescribed_dstar`             |
| `xblsys.f:BLKIN`  | `src/bl/station.rs::StationState::set_kinematic_variables` |
| `xblsys.f:BLVAR`  | `src/bl/station.rs::StationState::set_closure_variables`   |
| `xblsys.f:BLSYS`  | `src/bl/blsys.rs::assemble_interval_system`            |
| `xblsys.f:BLDIF`  | `src/bl/difference.rs::IntervalSystem::assemble_interval_equations` |
| `xblsys.f:TRCHEK2`| `src/bl/transition.rs::check_transition`               |
| `xsolve.f:BLSOLV` | `src/bl/blsolv.rs::solve_newton_system`                |
| `xpanel.f:PSILIN` | `src/solver/psilin.rs::panel_influence`                |
| `xpanel.f:GGCALC` | `src/solver/ggcalc.rs::build_inviscid_system`          |

Every translating function carries `#[doc(alias = "XFOIL NAME")]`, so `cargo doc` search by the
Fortran name finds it.

## Pages of the translation map

| Page | What it holds |
| --- | --- |
| [Variable and routine names](variable-and-routine-names.md) | XFOIL name to yFoil name, for every translated symbol, table by table |
| [Boundary-layer state layout](boundary-layer-state-layout.md) | How yFoil's station arrays mirror XFOIL's two sides with the wake on side 2 |
| [Architecture diagram](architecture-diagram.md) | XFOIL's solution method as a diagram: geometry, inviscid solution, boundary layer, coupling |
| [XFOIL's control flow and variables](xfoils-control-flow-and-variables.md) | Which routine sets which COMMON variables, in a typical run |
| [XFOIL's source files, routine by routine](xfoils-source-files-routine-by-routine.md) | A walk through the translated routines of `spline.f`, `xpanel.f`, `xbl.f`, `xblsys.f` and `xoper.f` |
| [XFOIL's COMMON blocks](xfoils-common-blocks.md) | `XFOIL.INC`, `XBL.INC` and `BLPAR.INC`, block by block |

The first two are what most readers of this tab need; the rest is background on XFOIL itself,
written during the translation, and is consulted rather than read.
