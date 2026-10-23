# 1.5 Formatted output precision — Overcome (instrumentation)

`PSAVE` is `G15.7`, `.pol` is `F9.4/F10.5`, `DUMP` is `F9.5/F10.6`, `CPWR` is `F11.5`: none can
support a comparison below ~1e-5, and XFOIL's stock debug formats `E24.16/E18.10/E12.4` carry 16, 10
and 4 significant figures. The instrumented reference build writes `ES24.16` (17 figures, the minimum
that round-trips an f64); the instrumentation is proven inert (byte-identical pristine vs instrumented
outputs, `scripts/xfoil-build.sh --verify`). No test compares against XFOIL's formatted files (CLAUDE.md
Rule 4). An old validation table that showed 4.9e-4 to 1.6e-3 RMS "paneling error" was measuring
7-digit PSAVE output plus the two NACA definitions of §6.1 (`docs/development/testing/geometry/README.md`).
