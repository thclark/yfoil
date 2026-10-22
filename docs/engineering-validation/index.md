---
icon: lucide/scale
---

# Engineering validation

This tab is for engineers deciding whether yFoil's results are fit for aerodynamic work. It
answers that the way an aerodynamicist would check any code: are the sections the ones the
designations name, do the results obey the physics, how do they compare with the code everyone
already trusts on the cases people actually run, and how much do they move when the inputs
change slightly. The results read as polars, transition locations and section shapes.

It does not ask whether yFoil is a faithful *translation* of XFOIL, statement by statement — that
is a developer's question, with its own comparison target (the instrumented XFOIL itself), and it
is the next tab, [Fidelity to XFOIL](../fidelity-to-xfoil/index.md). The two meet in one place:
this tab compares whole sweeps with XFOIL under normal conditions, as an engineer reads them; the
next holds yFoil to XFOIL one step at a time, which is what the tests check.

## Pages of Engineering validation, in reading order

| Page | Question it answers | Compared against |
| --- | --- | --- |
| [The validation set](the-validation-set.md) | Which sections and flow conditions are assessed, and why each | — |
| [Generated sections against naca456](generated-sections-against-naca456.md) | Is each generated section the one its designation names? | NASA's public-domain `naca456` ordinates |
| [Physical invariants](physical-invariants.md) | Do the results obey what the physics requires, whatever XFOIL does? | symmetry, thin-aerofoil theory, Blasius |
| [Agreement with XFOIL under normal conditions](agreement-with-xfoil-under-normal-conditions.md) | Does yFoil give XFOIL's polars across the sections and regimes people run? | the reference XFOIL, point by point |
| [Sensitivity of results to their inputs](sensitivity-of-results-to-their-inputs.md) | How far do XFOIL's results move when the inputs change slightly? | the reference itself, perturbed |

## How these pages are made

The generated pages (naca456, agreement with XFOIL) are written by the study that computed them
and can be regenerated at any time from tracked inputs; nothing on them is hand-fed, and no
reference data is committed under `docs/`. The command for each is on
[regenerating the reports](../contributing/regenerating-the-reports.md).

Figures are [ironlab](https://github.com/thclark/ironlab) figures, built in Rust by each study
from the numbers it computes ([figures](../contributing/conventions/figures.md)). They are
embedded live in these pages: pan, zoom and toggle lines from the legend. Each carries parameters
(case, conditions, outputs) and labels, so every study's figures can be loaded together into the
ironlab viewer and sorted, filtered and grouped. A study's run folder also holds a PDF of each.

**Not here yet.** Boundary-layer distributions — displacement and momentum thickness and skin
friction at selected angles, both codes compared — for the ±15° N = 160 polar are the remaining
validation stage and will land on the agreement page. An earlier version carried XFOIL's
formatted `DUMP`/`CPWR`/`PACC` files (F9.5, F11.5, F10.5 — about 1e-5) as committed reference
data, which cannot support the project's tolerances, and was removed on 2026-09-04.
