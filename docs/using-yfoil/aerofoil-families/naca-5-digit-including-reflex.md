# NACA 5-digit, including reflex

*Report 537* ([jacobs1935](../bibliography.md#jacobs1935)) put the maximum camber unusually far
forward to raise the maximum lift with little pitching moment. Designation `LPRtt`: the design
lift coefficient \(c_l = 0.15\,L\), the position of maximum camber \(p = P/20\), \(R = 0\) for the
standard line or 1 for the reflexed one, and \(t\) in percent. The thickness form is the
[4-digit one](naca-4-digit-and-4-digit-modified.md#naca-4-digit).

## The 3-digit mean line

The **3-digit mean line** is a cubic forward of \(x = r\) and a straight line aft:

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

\(k_1\) scales linearly with \(c_l/0.3\).

## The reflexed mean line

The **reflexed line** adds a second cubic aft, with
\(k_2/k_1 = \left(3(r - p)^2 - r^3\right)/(1 - r)^3\) chosen for zero pitching moment
(Report 537's expression; both NASA memoranda misprint it, as
[carmichael2001](../bibliography.md#carmichael2001) notes):

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

## The 5-digit family drawn

The 230 line (23012, 23015, 23018) is on a great many aircraft, from the Cessna singles to the
Beech King Air ([lednicer](../bibliography.md#lednicer)).

<ironlab-figure src="naca-5-digit-including-reflex/naca-5-digit.fig" format="fig" name="naca-5-digit" alt="NACA 5-digit: default NACA 23018 in black; design lift coefficient cl (p = 0.15, t = 0.18); maximum-camber position p (cl = 0.3, t = 0.18); reflex mean line, position p (cl = 0.3, t = 0.18); maximum thickness t (230 mean line)">NACA 5-digit: default NACA 23018 in black; design lift coefficient cl (p = 0.15, t = 0.18); maximum-camber position p (cl = 0.3, t = 0.18); reflex mean line, position p (cl = 0.3, t = 0.18); maximum thickness t (230 mean line)</ironlab-figure>

*NACA 5-digit: default NACA 23018 in black; design lift coefficient cl (p = 0.15, t = 0.18); maximum-camber position p (cl = 0.3, t = 0.18); reflex mean line, position p (cl = 0.3, t = 0.18); maximum thickness t (230 mean line).*

The generator is checked against naca456 in
[NACA 5-digit against naca456](../../engineering-validation/generated-sections-against-naca456.md#naca-5-digit-against-naca456).
