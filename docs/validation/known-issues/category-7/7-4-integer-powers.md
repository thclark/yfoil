# 7.4 Integer powers — Replicated (gfortran)

gfortran expands `X**2`, `X**3` inline as multiplications; yFoil writes `x*x*x` where bit-exactness
matters rather than `powi` (CLAUDE.md conventions). Transcendentals are not the host libm's in either
code: the reference is linked against CORE-MATH's correctly rounded functions and yFoil calls the same
ones (`cr-libm`), so their bits do not depend on the host. Before 2026-10-15 they were the host libm's:
Apple libSystem and glibc differ by 1 ULP on 0.1 % (`exp`, `ln`, `pow`) to 18 % (`tanh`) of inputs
(measured 2026-09-11), and XFOIL's own post-CL_max wanderings differed by O(1) between the two while its
converged points moved by ≤ 1.4e-11 (`docs/cross-platform-quirks.md`).
