# NACA 6- and 6A-series

The laminar-flow sections, whose thickness forms survive only as tables, and their revision with
straight-sided trailing edges.

## NACA 6-series

*Report 824* ([abbott1945](../bibliography.md#abbott1945)) and *Theory of Wing Sections*
([abbott1959](../bibliography.md#abbott1959)): the laminar-flow sections. Designation `6F-Ltt`
(`63-415`): the family digit \(F\) is the position of minimum pressure in tenths of chord, \(L\)
the design lift coefficient in tenths, \(t\) in percent. A subscript or bracketed digit after the
family (`64(1)-212`) states the low-drag range and does not change the geometry; yFoil ignores it.

### The 6-series thickness forms

**The thickness forms have no closed form.** They were derived by conformal mapping
([theodorsen1931](../bibliography.md#theodorsen1931),
[theodorsen1933](../bibliography.md#theodorsen1933)) for a prescribed pressure distribution, and
survive as 201-point tables of the mapping functions \(\varepsilon(\phi)\) and \(\psi(\phi)\)
for each of 63, 64, 65, 66, 67, 63A, 64A and 65A
([ladson1974](../bibliography.md#ladson1974), [ladson1996](../bibliography.md#ladson1996)). With
\(s(t)\) the scale factor that gives thickness ratio \(t\) (a quartic fit per family),

\[
z = e^{\,s\psi_0 + i\phi},\qquad
z' = z\,e^{\,s(\psi - \psi_0) - i\,s\varepsilon},\qquad
\zeta = z' + \frac{1}{z'},\qquad
x + iy = \frac{\zeta_0 - \zeta}{|\zeta_{200} - \zeta_0|}
\]

gives the surface from the leading edge (\(\phi = 0\)) to the trailing edge (\(\phi = \pi\)),
which closes exactly. Ordinates at other stations come from the arc-length cubic spline through
the 201 points. TM-4741 records that the original graphs of \(\varepsilon\) and \(\psi\) are
lost and that the program reproduces the published ordinates to within \(5 \times 10^{-5}\)
chord; the tables are the definition.

### The a mean line

The **a mean line** carries uniform load from the leading edge to \(x = a\), falling linearly to
zero at the trailing edge ([abbott1959](../bibliography.md#abbott1959), eq. 4.26):

\[
y_c = \frac{c_l}{2\pi(a + 1)}\left\{
\frac{1}{1 - a}\left[\tfrac{1}{2}(a - x)^2\ln|a - x| - \tfrac{1}{2}(1 - x)^2\ln(1 - x)
+ \tfrac{1}{4}(1 - x)^2 - \tfrac{1}{4}(a - x)^2\right] - x\ln x + g - hx\right\}
\]

\[
g = -\frac{1}{1 - a}\left[a^2\left(\tfrac{1}{2}\ln a - \tfrac{1}{4}\right) + \tfrac{1}{4}\right],\qquad
h = \frac{1}{1 - a}\left[\tfrac{1}{2}(1 - a)^2\ln(1 - a) - \tfrac{1}{4}(1 - a)^2\right] + g
\]

and for \(a = 1\), the default and the line the designations imply,
\(y_c = -\dfrac{c_l}{4\pi}\left[(1 - x)\ln(1 - x) + x\ln x\right]\). Another \(a\) is
`yfoil geometry naca 63-415 --a 0.5`.

The 63-415 was the standard wind-turbine section before dedicated families appeared, and the
Risø catalogue ([bertagnolio2001](../bibliography.md#bertagnolio2001)) compares measurements with
XFOIL and Navier–Stokes results for it.

<ironlab-figure src="naca-6-and-6a-series/naca-6-series.fig" format="fig" name="naca-6-series" alt="NACA 6-series: default NACA 63-415 in black; family: minimum-pressure position (cl = 0.4, t = 0.15); design lift coefficient cl (63, t = 0.15); extent of uniform loading a (63-415); maximum thickness t (63, cl = 0.4)">NACA 6-series: default NACA 63-415 in black; family: minimum-pressure position (cl = 0.4, t = 0.15); design lift coefficient cl (63, t = 0.15); extent of uniform loading a (63-415); maximum thickness t (63, cl = 0.4)</ironlab-figure>

*NACA 6-series: default NACA 63-415 in black; family: minimum-pressure position (cl = 0.4, t = 0.15); design lift coefficient cl (63, t = 0.15); extent of uniform loading a (63-415); maximum thickness t (63, cl = 0.4).*

## NACA 6A-series

*Report 903* ([loftin1948](../bibliography.md#loftin1948)) removed the 6-series' cusped trailing
edge: the 6A thickness forms are straight-sided from 0.8 chord, and the **6A mean line** is the
\(a = 0.8\) line scaled by 0.97948 and replaced by a straight line aft of 0.86 chord, so that the
cambered sections are straight-sided too:

\[
y_c = \begin{cases}
0.97948\;y_c^{(a = 0.8)}(x) & x \le 0.86 \\
-0.24521\,c_l\,(x - 1) & x > 0.86
\end{cases}
\]

Designation `6FALtt` (`64A010`, `63A415`). The 64A010 is the standard transonic and
unsteady-aerodynamics test section.

<ironlab-figure src="naca-6-and-6a-series/naca-6a-series.fig" format="fig" name="naca-6a-series" alt="NACA 6A-series: default NACA 64A010 in black; family (cl = 0, t = 0.10); design lift coefficient cl (64A, t = 0.10); maximum thickness t (64A, cl = 0)">NACA 6A-series: default NACA 64A010 in black; family (cl = 0, t = 0.10); design lift coefficient cl (64A, t = 0.10); maximum thickness t (64A, cl = 0)</ironlab-figure>

*NACA 6A-series: default NACA 64A010 in black; family (cl = 0, t = 0.10); design lift coefficient cl (64A, t = 0.10); maximum thickness t (64A, cl = 0).*

Both generators are checked against naca456 in
[NACA 6-series against naca456](../../engineering-validation/generated-sections-against-naca456.md#naca-6-series-against-naca456)
and the section after it; the tabulated thickness forms are where the two differ most, by
naca456's own root tolerance.
