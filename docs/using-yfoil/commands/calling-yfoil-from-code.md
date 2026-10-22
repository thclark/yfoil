---
icon: lucide/plug
---

# Calling yFoil from code

The commands are a thin layer over a Rust library, and the library is the only other way
to call yFoil today. **There are no Python, Julia or WebAssembly bindings yet**; they are
planned, and this page will describe them when they exist.

## Calling the Rust library

[Install](../installing-yfoil.md#as-a-rust-library) the crate, then generate a section, panel it
and solve a point. This is the repository's `examples/validate_viscous.rs`, shortened:

``` rust
use yfoil::geometry::Thickness;
use yfoil::geometry::{naca_4digit, panel_foil};
use yfoil::solver::analysis::{analyse, FlowConditions};

fn main() {
    let geometry = naca_4digit("0012", 160, Thickness::Perpendicular).expect("NACA 0012");
    let airfoil = panel_foil(&geometry);
    let spec = FlowConditions {
        re: Some(1.0e6),
        mach: 0.0,
        ncrit: 9.0,
        max_iterations: 20,
        ..FlowConditions::default()
    };
    let p = analyse(&airfoil, 4.0_f64.to_radians(), &spec);
    println!("CL {:.6}  CD {:.6}  converged {}", p.cl, p.cd, p.converged);
}
```

`analyse` is one operating point from scratch, the library's [`yfoil analyse`](analysing-one-operating-point.md).
A sequence of points that each start from the last, as XFOIL's OPER does, is a
`yfoil::solver::analysis::Session`: `Session::new(&airfoil, spec)`, then `alpha`, `cl` and
`init` in the order XFOIL's commands would come; `compute_polar` is the whole
[polar sweep](sweeping-a-polar.md). Angles are in radians in the library and in degrees on the
command line.

Inside `solver/` and `bl/` the library mirrors XFOIL's data model exactly, which is a
deliberate choice: it keeps every translated statement diffable against its Fortran original.
Where each XFOIL variable and routine went is the
[translation map](../../fidelity-to-xfoil/translation-map/index.md).
