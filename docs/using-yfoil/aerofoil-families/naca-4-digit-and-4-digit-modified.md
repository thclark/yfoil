# NACA 4-digit and 4-digit modified

The oldest family and its variant with a separately chosen nose and thickness position. Both use
the 2-digit mean line; the modified family replaces the thickness form.

## NACA 4-digit

*Report 460* ([jacobs1933](../bibliography.md#jacobs1933)) fitted a four-term polynomial to the
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

<ironlab-figure src="naca-4-digit-and-4-digit-modified/naca-4-digit.fig" format="fig" name="naca-4-digit" alt="NACA 4-digit: default NACA 0012 in black; maximum camber m (p = 0.4, t = 0.12); maximum-camber position p (m = 0.04, t = 0.12); maximum thickness t (m = 0)">NACA 4-digit: default NACA 0012 in black; maximum camber m (p = 0.4, t = 0.12); maximum-camber position p (m = 0.04, t = 0.12); maximum thickness t (m = 0)</ironlab-figure>

*NACA 4-digit: default NACA 0012 in black; maximum camber m (p = 0.4, t = 0.12); maximum-camber position p (m = 0.04, t = 0.12); maximum thickness t (m = 0).*

## NACA 4-digit modified

*Report 492* ([stack1934](../bibliography.md#stack1934)) varied the leading-edge radius and the
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
its slope there ([ladson1975](../bibliography.md#ladson1975), TM-4741 p. 6). The mean line is the
2-digit line. The `0012-63` is the 4-digit 0012 to a few \(10^{-4}\); the family separates the
leading-edge radius from the thickness, which is what it is for.

<ironlab-figure src="naca-4-digit-and-4-digit-modified/naca-4-digit-modified.fig" format="fig" name="naca-4-digit-modified" alt="NACA 4-digit modified: default NACA 0012-34 in black; leading-edge radius index I (M = 0.4, t = 0.12); maximum-thickness position M (I = 3, t = 0.12); maximum thickness t (I = 3, M = 0.4)">NACA 4-digit modified: default NACA 0012-34 in black; leading-edge radius index I (M = 0.4, t = 0.12); maximum-thickness position M (I = 3, t = 0.12); maximum thickness t (I = 3, M = 0.4)</ironlab-figure>

*NACA 4-digit modified: default NACA 0012-34 in black; leading-edge radius index I (M = 0.4, t = 0.12); maximum-thickness position M (I = 3, t = 0.12); maximum thickness t (I = 3, M = 0.4).*

Both generators are checked against naca456 in
[NACA 4-digit against naca456](../../engineering-validation/generated-sections-against-naca456.md#naca-4-digit-against-naca456)
and the section after it.
