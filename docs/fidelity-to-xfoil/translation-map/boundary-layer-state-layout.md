# Boundary-layer state layout

XFOIL uses COMMON blocks; yFoil passes state explicitly. Everywhere else in the crate that frees
the Rust to be idiomatic, but inside `src/solver/` and `src/bl/` the boundary-layer state mirrors
XFOIL's data model **exactly**, so that a translated statement corresponds one for one to its
Fortran line and a dumped XFOIL array can be loaded into yFoil without reindexing.

## Two sides, the wake on side 2

- There are **two sides**: side 1 is the upper surface and side 2 the lower, each marched from the
  stagnation point to the trailing edge.
- **The wake is appended to side 2**: `NBL(2) = IBLTE(2) + NW`, so side 2's stations run on past
  the trailing-edge station into the wake. There is no third "wake surface".
- Station arrays are **1-based with a dummy slot 0**, indexed `[side][i_station]`, so
  `theta[2][i]` is XFOIL's `THET(I,2)`.

## The maps between the frames

Four discretisations run at once — buffer geometry, panel nodes, BL stations and Newton rows (the
suffixes are in [naming](../../contributing/conventions/naming.md#index-systems)) — and the
boundary-layer state carries the maps between the last three explicitly, as XFOIL does:

| yFoil | XFOIL | Maps |
| --- | --- | --- |
| `i_node[side][i_station]` | `IPAN(IBL,IS)` | station → the panel node it sits on |
| `i_row[side][i_station]` | `ISYS(IBL,IS)` | station → its row of the global Newton system |
| `velocity_sign[side][i_station]` | `VTI(IBL,IS)` | ±1 between the panel's tangential velocity and the BL edge velocity |
| `i_te_station[side]` | `IBLTE(IS)` | the trailing-edge station of each side |
| `n_stations[side]` | `NBL(IS)` | the last station of each side, side 2 including the wake |
| `i_transition_station[side]` | `ITRAN(IS)` | the station of the transition interval |

The rest of the per-station state — `ue`, `theta`, `dstar`, `sqrtctau`, `mass_defect` and their
closures — is listed with its XFOIL name in Table 2 of
[variable and routine names](variable-and-routine-names.md);
XFOIL's own declaration of the same arrays is in
[XFOIL's COMMON blocks](xfoils-common-blocks.md).

## Short locals inside a translated body

Inside a translated routine, short locals are bound at the top from the named fields
(`let hk1 = station1.hk;`), so the formula lines below them stay diffable against the Fortran
while the fields keep their descriptive names.
