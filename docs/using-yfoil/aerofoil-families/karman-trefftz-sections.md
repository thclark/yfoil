# Kármán–Trefftz sections

Not a NACA family: the conformal map of a circle ([karman1918](../bibliography.md#karman1918)),
generalising Joukowski's transformation ([joukowski1910](../bibliography.md#joukowski1910)) to a
finite trailing-edge angle \(\tau\). A circle through \(\zeta = 1\) with centre
\((x_c, y_c)\) is mapped by

\[
z = n\,\frac{1 + w^n}{1 - w^n},\qquad w = \frac{\zeta - 1}{\zeta + 1},\qquad n = 2 - \frac{\tau}{\pi}
\]

(\(\tau = 0\) is the Joukowski map \(z = \zeta + 1/\zeta\) and its cusp). \(x_c < 0\) sets the
thickness, \(y_c > 0\) the camber. The image is normalised so the trailing edge (the image of
\(\zeta = 1\)) is at \(x = 1\) and the leading edge (minimum \(x\)) at \(x = 0\).

``` sh
yfoil geometry karman-trefftz --x-centre -0.1 --y-centre 0.05 --te-angle 10 -n 160 -o kt.json
```

The section is analytic: its potential-flow solution is known in closed form, so it checks a
panel method against an exact answer with no boundary layer in the way, and its sharp trailing
edge with a finite angle exercises the solver's `SHARP` branch. It is not a wing section anyone
builds. naca456 does not generate it, so it has no comparison there.

<ironlab-figure src="karman-trefftz-sections/karman-trefftz.fig" format="fig" name="karman-trefftz" alt="Kármán–Trefftz: default Karman-Trefftz xc=-0.1 yc=0.05 tau=10 in black; circle centre x (yc = 0.05, τ = 10°); circle centre y (xc = −0.1, τ = 10°); trailing-edge angle τ (xc = −0.1, yc = 0.05)">Kármán–Trefftz: default Karman-Trefftz xc=-0.1 yc=0.05 tau=10 in black; circle centre x (yc = 0.05, τ = 10°); circle centre y (xc = −0.1, τ = 10°); trailing-edge angle τ (xc = −0.1, yc = 0.05)</ironlab-figure>

*Kármán–Trefftz: default Karman-Trefftz xc=-0.1 yc=0.05 tau=10 in black; circle centre x (yc = 0.05, τ = 10°); circle centre y (xc = −0.1, τ = 10°); trailing-edge angle τ (xc = −0.1, yc = 0.05).*
