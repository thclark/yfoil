---
icon: lucide/scale
---

# Against naca456

yFoil's NACA generators are checked against the public-domain NASA/PDAS `naca456` program
(Carmichael 2001, the revision of the NASA Langley ordinate programs of TM X-3284, TM X-3069 and
TM-4741), run by `scripts/naca456-fixtures.sh` through a driver on its own modules at 17 significant
figures. The comparison is the evidence that each generator produces the published section, which
matters because every case on these pages starts from yFoil's panels (XFOIL loads them rather than
generating its own).

## How the difference is measured

Differences are `|Δ| / max(|yFoil|, |naca456|, 1)`, the worst over naca456's 98 "very fine"
stations, for the thickness form y_t and its slope, the mean line y_c and its slope, and the surface.

## Results

--8<-- "validation/aerofoil-series/tables/naca456.md"

## What the differences mean

The closed-form families agree to round-off. The 6-series mean lines differ by naca456's 10-digit π
(1.1e-10 relative). The 6-series thickness forms differ by naca456's own root tolerance: it finds the
station on the arc-length spline with Brent's method at 1e-6 and reports the ordinate at the station
reached, so its ordinates are off by up to 1e-6 × |dy_t/dx| (largest at the nose); yFoil inverts to
round-off. The tolerances in `tests/application/naca456.rs` are derived from exactly that
(`tests/common/utilities/tolerances.rs`). naca456 reports the aft slope of the 4-digit modified form
as dy/d(1 − x); the sign is restored before comparing.
