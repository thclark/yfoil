# NACA 64A010, TYPE 2

Shows: [§7.13 Whether a sweep converges, and what it converges to, can depend on the panel count](../category-7/7-13-whether-a-sweep-converges-and-what-it-converges-to-can-depend-on-the-panel-count.md); [§4 Silent recovery](../category-4/index.md).

--8<-- "validation/known-issues/tables/cases-naca-64a010-type-2.md"

!!! note "Stub — to be written"

    What this run shows, and where and why yFoil and XFOIL part.

## N = 60

`naca64a010_n60_polar30_re1e6_type2` — Known issues 7.13 at N = 60: NACA 64A010, viscous TYPE 2, M 0 (MRCL 100xRe clamp at CL ~ 0; 5 of 14 points converge).

<ironlab-figure src="../figures/polars-naca64a010_n60_polar30_re1e6_type2.fig" format="fig" name="polars-naca64a010_n60_polar30_re1e6_type2" alt="naca64a010_n60_polar30_re1e6_type2: CL, CD, the last RMSBL, iterations, transition location and position in its station interval against alpha, XFOIL and yFoil">naca64a010_n60_polar30_re1e6_type2: CL, CD, the last RMSBL, iterations, transition location and position in its station interval against alpha, XFOIL and yFoil</ironlab-figure>

## N = 160

`naca64a010_n160_polar30_re1e6_type2_iter100` — Known issues 7.13 at N = 160: NACA 64A010, viscous TYPE 2, M 0 (11 of 13 converge; non-finite at 12 deg and hangs, 4).

<ironlab-figure src="../figures/polars-naca64a010_n160_polar30_re1e6_type2_iter100.fig" format="fig" name="polars-naca64a010_n160_polar30_re1e6_type2_iter100" alt="naca64a010_n160_polar30_re1e6_type2_iter100: CL, CD, the last RMSBL, iterations, transition location and position in its station interval against alpha, XFOIL and yFoil">naca64a010_n160_polar30_re1e6_type2_iter100: CL, CD, the last RMSBL, iterations, transition location and position in its station interval against alpha, XFOIL and yFoil</ironlab-figure>

## N = 240

`naca64a010_n240_polar30_re1e6_type2_iter200` — Known issues 7.13 at N = 240: NACA 64A010, viscous TYPE 2, M 0 (no point converges).

<ironlab-figure src="../figures/polars-naca64a010_n240_polar30_re1e6_type2_iter200.fig" format="fig" name="polars-naca64a010_n240_polar30_re1e6_type2_iter200" alt="naca64a010_n240_polar30_re1e6_type2_iter200: CL, CD, the last RMSBL, iterations, transition location and position in its station interval against alpha, XFOIL and yFoil">naca64a010_n240_polar30_re1e6_type2_iter200: CL, CD, the last RMSBL, iterations, transition location and position in its station interval against alpha, XFOIL and yFoil</ironlab-figure>
