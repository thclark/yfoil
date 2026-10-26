---
icon: lucide/shapes
---

# Generator families

One figure per generator family yFoil offers. Each panel varies one parameter: the family's member in the
validation set is drawn in black on top, and the members that vary the parameter sit behind it, coloured
by value (YlGnBu, lightest = smallest). Every panel is on the same linked axes (x from −0.2 to 1.2, y from
−0.5 to 0.5, at equal scale), so shapes can be compared across panels and across families. Pan and zoom
in one panel moves the others with it.

The formulas and history of each family are in the
[aerofoil series guide](../../user-guide/reference/aerofoil-series.md); the references are in the
[bibliography](../../user-guide/reference/bibliography.md).

## NACA 4-digit

Parameters: m (maximum camber), p (its position), t (maximum thickness).

The 1933 family of Report 460: a four-term polynomial thickness form on a two-parabola mean line, chosen to fit the best sections of the day (Clark Y, Göttingen 398). It remains the most used section family in general aviation and in codes' self-tests, and XFOIL's `NACA` command generates it.

<ironlab-figure src="figures/naca-4-digit.fig" format="fig" name="naca-4-digit" alt="NACA 4-digit: default NACA 0012 in black; maximum camber m (p = 0.4, t = 0.12); position of maximum camber p (m = 0.04, t = 0.12); maximum thickness t (m = 0)">NACA 4-digit: default NACA 0012 in black; maximum camber m (p = 0.4, t = 0.12); position of maximum camber p (m = 0.04, t = 0.12); maximum thickness t (m = 0)</ironlab-figure>

*NACA 4-digit: default NACA 0012 in black; maximum camber m (p = 0.4, t = 0.12); position of maximum camber p (m = 0.04, t = 0.12); maximum thickness t (m = 0).*

**In the validation set:** NACA 0012 is the symmetric baseline (Rule 6's XFOIL-independent invariants: C_L = C_M = 0 at α = 0, mirror symmetry at −α) and the CI reference case; NACA 4412 is the cambered 4-digit member and a real wing section (Lednicer lists the Champion 7EC Traveler).

References: [jacobs1933](../../user-guide/reference/bibliography.md#jacobs1933), [abbott1959](../../user-guide/reference/bibliography.md#abbott1959).

## NACA 4-digit modified

Parameters: m, p (the 2-digit mean line), t, I (leading-edge radius index), M (position of maximum thickness).

Report 492 (1934) varied the leading-edge radius and the position of maximum thickness of the 4-digit form independently, at speeds up to and beyond the compressibility burble: index I = 6 is the 4-digit radius, 0 a sharp nose; M is in tenths of chord. The aft part is a cubic in (1 − x) with a fixed trailing-edge thickness.

<ironlab-figure src="figures/naca-4-digit-modified.fig" format="fig" name="naca-4-digit-modified" alt="NACA 4-digit modified: default NACA 0012-34 in black; leading-edge radius index I (M = 0.4, t = 0.12); position of maximum thickness M (I = 3, t = 0.12); maximum thickness t (I = 3, M = 0.4)">NACA 4-digit modified: default NACA 0012-34 in black; leading-edge radius index I (M = 0.4, t = 0.12); position of maximum thickness M (I = 3, t = 0.12); maximum thickness t (I = 3, M = 0.4)</ironlab-figure>

*NACA 4-digit modified: default NACA 0012-34 in black; leading-edge radius index I (M = 0.4, t = 0.12); position of maximum thickness M (I = 3, t = 0.12); maximum thickness t (I = 3, M = 0.4).*

**In the validation set:** NACA 0012-34 has a sharper nose (I = 3) and its maximum thickness further aft than the 0012, with the same t: the pair separates the effect of leading-edge radius on transition and stagnation-point behaviour from that of thickness.

References: [stack1934](../../user-guide/reference/bibliography.md#stack1934), [abbott1959](../../user-guide/reference/bibliography.md#abbott1959), [ladson1975](../../user-guide/reference/bibliography.md#ladson1975).

## NACA 5-digit

Parameters: cl (design lift coefficient), p (position of maximum camber), reflex, t.

Report 537 (1935) put the maximum camber unusually far forward — a cubic mean line joined to a straight (or, in the reflexed lines, a second cubic) aft part — to raise the maximum lift with little pitching moment. The 230 line (cl = 0.3, p = 0.15) is by far the most used.

<ironlab-figure src="figures/naca-5-digit.fig" format="fig" name="naca-5-digit" alt="NACA 5-digit: default NACA 23018 in black; design lift coefficient cl (p = 0.15, t = 0.18); position of maximum camber p (cl = 0.3, t = 0.18); reflex mean line, position p (cl = 0.3, t = 0.18); maximum thickness t (230 mean line)">NACA 5-digit: default NACA 23018 in black; design lift coefficient cl (p = 0.15, t = 0.18); position of maximum camber p (cl = 0.3, t = 0.18); reflex mean line, position p (cl = 0.3, t = 0.18); maximum thickness t (230 mean line)</ironlab-figure>

*NACA 5-digit: default NACA 23018 in black; design lift coefficient cl (p = 0.15, t = 0.18); position of maximum camber p (cl = 0.3, t = 0.18); reflex mean line, position p (cl = 0.3, t = 0.18); maximum thickness t (230 mean line).*

**In the validation set:** NACA 23018 is a real aircraft section (Lednicer lists the Beech 90 King Air root; the 23012 is its tip) and, with the camber concentrated at the nose, loads the leading edge heavily: transition and separation are exercised at the nose rather than mid-chord.

References: [jacobs1935](../../user-guide/reference/bibliography.md#jacobs1935), [abbott1959](../../user-guide/reference/bibliography.md#abbott1959).

## NACA 16-series

Parameters: cl (design lift coefficient of the a = 1 mean line), t.

Stack's 1943 sections (Report 763) were designed to delay the compressibility burble: the 4-digit modified form with I = 4 and maximum thickness at mid-chord (a 16-0tt is a 00tt-45), on the uniform-load a = 1 mean line. They were developed for propellers, and 24 of them were tested at Mach 0.3–0.8 (TN 1546).

<ironlab-figure src="figures/naca-16.fig" format="fig" name="naca-16" alt="NACA 16-series: default NACA 16-212 in black; design lift coefficient cl (t = 0.12); maximum thickness t (cl = 0.2)">NACA 16-series: default NACA 16-212 in black; design lift coefficient cl (t = 0.12); maximum thickness t (cl = 0.2)</ironlab-figure>

*NACA 16-series: default NACA 16-212 in black; design lift coefficient cl (t = 0.12); maximum thickness t (cl = 0.2).*

**In the validation set:** NACA 16-212 is the compressible-flow member of the set: a section designed and measured at subcritical Mach numbers, where XFOIL's Kármán–Tsien correction applies.

References: [stack1943](../../user-guide/reference/bibliography.md#stack1943), [lindsey1948](../../user-guide/reference/bibliography.md#lindsey1948), [ladson1975](../../user-guide/reference/bibliography.md#ladson1975).

## NACA 6-series

Parameters: family (63…67: position of minimum pressure), cl (design lift coefficient), a (extent of uniform loading), t.

The laminar-flow sections of Report 824 (1945): thickness forms derived by conformal mapping for a prescribed pressure distribution — the family digit is the position of minimum pressure in tenths of chord — on the a mean line, which carries uniform load to x = a and falls linearly to zero at the trailing edge. There is no closed form: the sections are defined by the tabulated ε and ψ mapping functions of TM-4741.

<ironlab-figure src="figures/naca-6.fig" format="fig" name="naca-6" alt="NACA 6-series: default NACA 63-415 in black; family: position of minimum pressure (cl = 0.4, t = 0.15); design lift coefficient cl (63, t = 0.15); extent of uniform loading a (63-415); maximum thickness t (63, cl = 0.4)">NACA 6-series: default NACA 63-415 in black; family: position of minimum pressure (cl = 0.4, t = 0.15); design lift coefficient cl (63, t = 0.15); extent of uniform loading a (63-415); maximum thickness t (63, cl = 0.4)</ironlab-figure>

*NACA 6-series: default NACA 63-415 in black; family: position of minimum pressure (cl = 0.4, t = 0.15); design lift coefficient cl (63, t = 0.15); extent of uniform loading a (63-415); maximum thickness t (63, cl = 0.4).*

**In the validation set:** NACA 63-415 was the standard wind-turbine section before dedicated families appeared; the Risø catalogue (Bertagnolio et al. 2001) compares measurements with both Navier–Stokes and XFOIL results for it. Its trailing edge closes, so it runs XFOIL's sharp-TE branch, and the 6-series is the family whose generator is *not* a formula.

References: [abbott1945](../../user-guide/reference/bibliography.md#abbott1945), [abbott1959](../../user-guide/reference/bibliography.md#abbott1959), [ladson1974](../../user-guide/reference/bibliography.md#ladson1974), [ladson1996](../../user-guide/reference/bibliography.md#ladson1996), [carmichael2001](../../user-guide/reference/bibliography.md#carmichael2001), [bertagnolio2001](../../user-guide/reference/bibliography.md#bertagnolio2001).

## NACA 6A-series

Parameters: family (63A, 64A, 65A), cl (design lift coefficient of the 6A mean line), t.

Loftin's 1948 revision (Report 903) of the 6-series removed the cusped trailing edge: the 6A thickness forms are straight-sided from 0.8 chord, and the modified a = 0.8 mean line is straight aft of 0.86 chord, so the cambered sections are too. The trailing-edge angle is finite and the sections are easier to build.

<ironlab-figure src="figures/naca-6a.fig" format="fig" name="naca-6a" alt="NACA 6A-series: default NACA 64A010 in black; family (cl = 0, t = 0.10); design lift coefficient cl (64A, t = 0.10); maximum thickness t (64A, cl = 0)">NACA 6A-series: default NACA 64A010 in black; family (cl = 0, t = 0.10); design lift coefficient cl (64A, t = 0.10); maximum thickness t (64A, cl = 0)</ironlab-figure>

*NACA 6A-series: default NACA 64A010 in black; family (cl = 0, t = 0.10); design lift coefficient cl (64A, t = 0.10); maximum thickness t (64A, cl = 0).*

**In the validation set:** NACA 64A010 is the classic transonic and unsteady-aerodynamics test section (the AGARD standard cases): symmetric, sharp trailing edge with a finite angle, straight-sided aft.

References: [loftin1948](../../user-guide/reference/bibliography.md#loftin1948), [ladson1974](../../user-guide/reference/bibliography.md#ladson1974), [ladson1996](../../user-guide/reference/bibliography.md#ladson1996).

## Kármán–Trefftz

Parameters: x_centre (thickness), y_centre (camber), τ (trailing-edge angle; 0 is the Joukowski cusp).

The conformal map of a circle (von Kármán and Trefftz 1918), generalising Joukowski's 1910 transformation to a finite trailing-edge angle. The potential-flow solution is known in closed form, so the section checks a panel method against an exact answer without a boundary layer in the way. It is not a wing section anyone builds.

<ironlab-figure src="figures/karman-trefftz.fig" format="fig" name="karman-trefftz" alt="Kármán–Trefftz: default Karman-Trefftz xc=-0.1 yc=0.05 tau=10 in black; circle centre x (yc = 0.05, τ = 10°); circle centre y (xc = −0.1, τ = 10°); trailing-edge angle τ (xc = −0.1, yc = 0.05)">Kármán–Trefftz: default Karman-Trefftz xc=-0.1 yc=0.05 tau=10 in black; circle centre x (yc = 0.05, τ = 10°); circle centre y (xc = −0.1, τ = 10°); trailing-edge angle τ (xc = −0.1, yc = 0.05)</ironlab-figure>

*Kármán–Trefftz: default Karman-Trefftz xc=-0.1 yc=0.05 tau=10 in black; circle centre x (yc = 0.05, τ = 10°); circle centre y (xc = −0.1, τ = 10°); trailing-edge angle τ (xc = −0.1, yc = 0.05).*

**In the validation set:** The analytic member: its exact inviscid pressure distribution is an XFOIL-independent reference for the panel solver, and its sharp trailing edge with a finite angle exercises the SHARP branch on a smooth, non-NACA shape.

References: [karman1918](../../user-guide/reference/bibliography.md#karman1918), [joukowski1910](../../user-guide/reference/bibliography.md#joukowski1910).
