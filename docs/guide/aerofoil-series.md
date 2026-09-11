---
icon: lucide/layers
---

# Aerofoil series

yFoil generates every NACA family that has a public, reproducible definition, and one analytic
section. A NACA section is a **basic thickness form** \(y_t(x)\) laid **perpendicular** to a
**mean line** \(y_c(x)\) ([abbott1959](../references.md#abbott1959), §6.2):

\[
x_u = x - y_t \sin\theta,\quad y_u = y_c + y_t \cos\theta,\qquad
x_l = x + y_t \sin\theta,\quad y_l = y_c - y_t \cos\theta,\qquad
\theta = \arctan\frac{\mathrm{d}y_c}{\mathrm{d}x}
\]

The formulas below are those of the NASA ordinate program
([ladson1996](../references.md#ladson1996)) and its public-domain revision `naca456`
([carmichael2001](../references.md#carmichael2001)), against which yFoil's generators are gated
(see [validation](../validation/aerofoil-series/README.md)). Every generated file carries a
[provenance record](#provenance) naming the family, the parameters and these references.

| Series | Designation | Thickness form | Mean line | Command |
| --- | --- | --- | --- | --- |
| [4-digit](#naca-4-digit) | `2412` | 4-digit | 2-digit | `yfoil geometry naca 2412` |
| [4-digit modified](#naca-4-digit-modified) | `0012-34`, `2412-63` | 4-digit modified | 2-digit | `yfoil geometry naca 0012-34` |
| [5-digit](#naca-5-digit) | `23012`, `23112` | 4-digit | 3-digit, reflex | `yfoil geometry naca 23012` |
| [16-series](#naca-16-series) | `16-212` | 4-digit modified (I = 4, M = 0.5) | 6-series, a = 1 | `yfoil geometry naca 16-212` |
| [6-series](#naca-6-series) | `63-415`, `64-010` | 63…67 (tabulated) | 6-series a | `yfoil geometry naca 63-415 [--a 0.5]` |
| [6A-series](#naca-6a-series) | `64A010`, `63A415` | 63A…65A (tabulated) | 6A | `yfoil geometry naca 64A010` |
| [Kármán–Trefftz](#karman-trefftz) | — | conformal map | — | `yfoil geometry karman-trefftz` |

All of them are sampled at the same cosine stations in \(x\) (\(n/2\) per surface, a node
exactly at the trailing edge, the leading edge straddled), so sections of different families are
compared on the same node distribution.

## NACA 4-digit

*Report 460* ([jacobs1933](../references.md#jacobs1933)) fitted a four-term polynomial to the
best sections of the day. Designation `mptt`: maximum camber \(m\) in percent, its position
\(p\) in tenths, maximum thickness \(t\) in percent of chord.

\[
y_t = 5t\left(0.2969\sqrt{x} - 0.1260\,x - 0.3516\,x^2 + 0.2843\,x^3 - 0.1015\,x^4\right)
\]

\[
y_c = \begin{cases}
\dfrac{m}{p^2}\left(2px - x^2\right) & x \le p \\[2ex]
\dfrac{m}{(1-p)^2}\left((1 - 2p) + 2px - x^2\right) & x > p
\end{cases}
\]

The \(-0.1015\) coefficient leaves a trailing-edge half-thickness of \(0.0105\,t\) (the section
is blunt; XFOIL's own `NACA` command uses the same coefficient); \(-0.1036\) would close it. The
family is the general-aviation workhorse and every code's self-test; XFOIL's canonical cases are
the 0012 and 4412.

## NACA 4-digit modified

*Report 492* ([stack1934](../references.md#stack1934)) varied the leading-edge radius and the
position of maximum thickness independently, at speeds up to and beyond the compressibility
burble. Designation `mptt-IM`: the 4-digit digits, then the **leading-edge radius index**
\(I\) (6 is the plain 4-digit radius, 0 a sharp nose, 9 three times the radius) and the
**position of maximum thickness** \(M\) in tenths of chord (2 to 6).

\[
y_t = \begin{cases}
a_0\sqrt{x} + a_1 x + a_2 x^2 + a_3 x^3 & x < M \\
d_0 + d_1(1-x) + d_2(1-x)^2 + d_3(1-x)^3 & x \ge M
\end{cases}
\qquad
r_{le} = 1.1019\left(\frac{tI}{6}\right)^2
\]

with \(d_0 = 0.002\) on the \(t = 0.2\) basis (the trailing-edge half-thickness, \(0.01\,t\)),
\(d_1\) the trailing-edge slope tabulated against \(M\), and the remaining coefficients fixed by
the leading-edge radius, the maximum thickness \(t\) at \(x = M\) and continuity of \(y_t\) and
its slope there ([ladson1975](../references.md#ladson1975), TM-4741 p. 6). The mean line is the
2-digit line. The `0012-63` is the 4-digit 0012 to a few \(10^{-4}\); the family separates the
leading-edge radius from the thickness, which is what it is for.

## NACA 5-digit

*Report 537* ([jacobs1935](../references.md#jacobs1935)) put the maximum camber unusually far
forward to raise the maximum lift with little pitching moment. Designation `LPRtt`: the design
lift coefficient \(c_l = 0.15\,L\), the position of maximum camber \(p = P/20\), \(R = 0\) for the
standard line or 1 for the reflexed one, and \(t\) in percent. The thickness form is the 4-digit
one. The **3-digit mean line** is a cubic forward of \(x = r\) and a straight line aft:

\[
y_c = \frac{k_1}{6}\begin{cases}
x^3 - 3rx^2 + r^2(3 - r)\,x & x < r \\
r^3(1 - x) & x \ge r
\end{cases}
\]

| \(p\) | 0.05 | 0.10 | 0.15 | 0.20 | 0.25 |
| --- | --- | --- | --- | --- | --- |
| \(r\) | 0.0580 | 0.1260 | 0.2025 | 0.2900 | 0.3910 |
| \(k_1\) (for \(c_l = 0.3\)) | 361.4 | 51.64 | 15.957 | 6.643 | 3.230 |

\(k_1\) scales linearly with \(c_l/0.3\). The **reflexed line** adds a second cubic aft, with
\(k_2/k_1 = \left(3(r - p)^2 - r^3\right)/(1 - r)^3\) chosen for zero pitching moment
(Report 537's expression; both NASA memoranda misprint it, as
[carmichael2001](../references.md#carmichael2001) notes):

\[
y_c = \frac{k_1}{6}\begin{cases}
(x - r)^3 - \dfrac{k_2}{k_1}(1 - r)^3 x - r^3 x + r^3 & x < r \\[2ex]
\dfrac{k_2}{k_1}(x - r)^3 - \dfrac{k_2}{k_1}(1 - r)^3 x - r^3 x + r^3 & x \ge r
\end{cases}
\]

| \(p\) | 0.10 | 0.15 | 0.20 | 0.25 |
| --- | --- | --- | --- | --- |
| \(r\) | 0.130 | 0.217 | 0.318 | 0.441 |
| \(k_1\) (for \(c_l = 0.3\)) | 51.99 | 15.793 | 6.520 | 3.191 |

The 230 line (23012, 23015, 23018) is on a great many aircraft, from the Cessna singles to the
Beech King Air ([lednicer](../references.md#lednicer)).

## NACA 16-series

*Report 763* ([stack1943](../references.md#stack1943)): sections designed to delay the
compressibility burble, developed for propellers. Designation `16-Ltt`: the design lift
coefficient \(c_l = L/10\) of a uniform-load (\(a = 1\)) mean line and \(t\) in percent. The
thickness form is the 4-digit modified form with \(I = 4\) and \(M = 0.5\), so a `16-012` is a
`0012-45`. Twenty-four of them were measured at Mach 0.3 to 0.8
([lindsey1948](../references.md#lindsey1948)), which makes the family the natural check of a
subcritical compressibility correction such as XFOIL's Kármán–Tsien.

## NACA 6-series

*Report 824* ([abbott1945](../references.md#abbott1945)) and *Theory of Wing Sections*
([abbott1959](../references.md#abbott1959)): the laminar-flow sections. Designation `6F-Ltt`
(`63-415`): the family digit \(F\) is the position of minimum pressure in tenths of chord, \(L\)
the design lift coefficient in tenths, \(t\) in percent. A subscript or bracketed digit after the
family (`64(1)-212`) states the low-drag range and does not change the geometry; yFoil ignores it.

**The thickness forms have no closed form.** They were derived by conformal mapping
([theodorsen1931](../references.md#theodorsen1931),
[theodorsen1933](../references.md#theodorsen1933)) for a prescribed pressure distribution, and
survive as 201-point tables of the mapping functions \(\varepsilon(\phi)\) and \(\psi(\phi)\)
for each of 63, 64, 65, 66, 67, 63A, 64A and 65A
([ladson1974](../references.md#ladson1974), [ladson1996](../references.md#ladson1996)). With
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

The **a mean line** carries uniform load from the leading edge to \(x = a\), falling linearly to
zero at the trailing edge ([abbott1959](../references.md#abbott1959), eq. 4.26):

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
Risø catalogue ([bertagnolio2001](../references.md#bertagnolio2001)) compares measurements with
XFOIL and Navier–Stokes results for it.

## NACA 6A-series

*Report 903* ([loftin1948](../references.md#loftin1948)) removed the 6-series' cusped trailing
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

## Kármán–Trefftz { #karman-trefftz }

Not a NACA family: the conformal map of a circle ([karman1918](../references.md#karman1918)),
generalising Joukowski's transformation ([joukowski1910](../references.md#joukowski1910)) to a
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
edge with a finite angle exercises the solver's `SHARP` branch.

## Trailing edges

The 4-digit families keep a finite trailing-edge thickness; the 6-series, 6A-series and
Kármán–Trefftz sections close exactly. XFOIL does not blunt a sharp trailing edge anywhere:
`PANGEN` keeps the buffer's end points, `TECALC` detects sharpness geometrically
(`SHARP = DSTE < 10⁻⁴ CHORD`) and the panel method, the wake and the boundary layer switch
branches on it, all of which yFoil translates and gates. yFoil therefore generates each section
as defined. Two adjustments are available, both recorded in the provenance record:

- `--sharp` moves the two trailing-edge nodes to their midpoint (a closed 4-digit section);
- `--te-gap GAP [--te-blend F]` is XFOIL's `TGAP`: the surfaces are moved apart along the
  existing gap's direction (or, for a closed edge, the mean end tangent) by
  \(\tfrac{1}{2}\Delta\,(x/c)\,e^{-(1 - x/c)(1/F - 1)}\) each, blended over the fraction \(F\)
  of the chord from the trailing edge. On a closed edge XFOIL delivers \(\cos(\tau/2)\) times the
  requested gap ([known issues §6.5](../xfoil-known-issues.md)); yFoil reproduces that.

## Provenance

Every generated geometry file carries a `generator` object: the series, the designation, the
thickness form and mean line with their parameter values, how the thickness was applied, whether
the trailing edge is sharp, any trailing-edge adjustment, the yFoil version, and the keys of the
references that define the family.

``` json
"generator": {
  "yfoil": "0.1.0",
  "series": "naca_6",
  "designation": "NACA 63-415",
  "thickness_form": { "family": "63", "t": 0.15 },
  "mean_line": { "family": "6", "cl": 0.4, "a": 1.0 },
  "thickness_applied": "perpendicular",
  "sharp_te": true,
  "references": ["abbott1945", "abbott1959", "ladson1974", "ladson1996", "carmichael2001"]
}
```

The record survives `yfoil geometry repanel` and `--te-gap` (they change the nodes, not the
section) and is absent from files converted from `.dat`.
