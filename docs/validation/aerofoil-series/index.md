# Aerofoil series

--8<-- "validation/aerofoil-series/tables/run.md"

## What this study is for

Every case in the validation runs on a section that yFoil generates and XFOIL loads, so both codes
see exactly the same panels. This study answers two questions about those sections: which ones the
cases use and why, and whether yFoil's generators produce the sections their designations define.

## The sections in the validation set

One section per generator family (two from the 4-digit series). Each was chosen for something the
others do not exercise, and each is public, mathematically defined and reproducible from its
designation.

--8<-- "validation/aerofoil-series/tables/selection.md"

## Main results

!!! note "Stub — to be written"

    The headline: every closed-form NACA generator matches naca456 to round-off, and the 6-series
    differs only by naca456's own root-finding tolerance.

## In this section

| Page | What it shows |
|---|---|
| [Generator families](generator-families.md) | Every family yFoil can generate, one figure per family and one panel per parameter, with why its member is in the validation set |
| [Against naca456](against-naca456.md) | The NACA generators compared with the NASA/PDAS reference program, and what the remaining differences mean |
