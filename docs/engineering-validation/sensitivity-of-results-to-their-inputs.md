# Sensitivity of results to their inputs

An engineer comparing two codes needs to know how much the reference itself moves when nothing
that matters has changed: if a node moved by a thousandth of a micron changes the answer, no
comparison can be sharper than that change. This page says how XFOIL's results respond to small
changes of three inputs — the node coordinates, the panel count and the alpha step — and where
each kind of evidence is.

## The input-sensitivity study

`scripts/input-sensitivity` runs the instrumented double-precision reference on one aerofoil
(NACA 0012 by default) over an upward alpha sweep, 0° to 25°, once for the base case and once per
perturbation level of three input families:

- **node coordinates**, every coordinate moved by 1 ULP up to 1e-7 chord;
- **panel count**, 160 down to 96;
- **alpha step**, 0.5° plus 1 ULP up to 1e-5°;

and records seven per-alpha quantities per level, drawn with the base case in front and the
perturbations behind it coloured by size. The node-perturbation method is in
`scripts/input-sensitivity/src/perturb.rs`.

``` sh
cargo run --release -p input-sensitivity
```

!!! warning "No published result yet"

    The study writes its figures and tables to a dated run folder under
    `scripts/input-sensitivity/runs/` (gitignored), not to this site. Its last run predates the
    move of both codes to correctly rounded maths (2026-10-15), so it measured a reference whose
    results also depended on the host's maths library. This page will carry its figures once the
    study is rerun and given a `--docs` mode.

## What is published

- **Every run's own 1-ULP twins.** Each case on
  [agreement with XFOIL under normal conditions](agreement-with-xfoil-under-normal-conditions.md)
  is run five more times with every panel coordinate jogged by one ULP; the page marks every
  point at which XFOIL was ill-conditioned — where its own twins part from it — so a difference
  from yFoil there is read as the reference's sensitivity, not as an error.
- **The panel count.** That a sweep's convergence, and what it converges to, can depend
  qualitatively on the panel count is measured at 60, 160 and 240 panels in
  [known issues §7.13](../fidelity-to-xfoil/xfoils-known-issues/the-issue-register.md), and shown
  in [runs that show the issues](../fidelity-to-xfoil/xfoils-known-issues/runs-that-show-the-issues.md).
- **The noise floor** from which the tests' tolerances are derived, stage by stage and variable by
  variable, is a developer's measurement:
  [the reference's noise floor](../fidelity-to-xfoil/tolerances/the-reference-noise-floor.md).
