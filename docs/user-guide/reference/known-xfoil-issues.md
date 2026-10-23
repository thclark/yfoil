---
icon: lucide/triangle-alert
---

# Known XFOIL issues

XFOIL is the reference (CLAUDE.md Rule 1) and its algorithm is reproduced including its own bugs and
quirks (Rule 2). This page is the register of everything found while translating it that a reader
comparing the two codes, or reading XFOIL's output, needs to know about. Each entry says what XFOIL
does, cites both sources, and states yFoil's handling:

- **Replicated** — yFoil does the same thing, and the behaviour is fixture-tested where a test exists.
- **Overcome** — yFoil deliberately differs. Only one such divergence affects results (the NACA
  generator); the others are output-side or instrumentation-side and leave the solver untouched.
- **Out of scope** — code XFOIL itself never executes on the analysis path, or design/plotting
  features yFoil does not translate. Listed so nobody looks for them.

Line numbers refer to `xfoil/third-party/xfoil-6.99/src/` and to yFoil at the commit that last
touched this page. Every "unreachable" claim below is also carried, with its reason, in
`xtask/fixtures-config/coverage.toml` and measured by `cargo xtask coverage`
(`docs/development/testing/coverage.md`).

## Categories

| | Category | Description |
|---|---|---|
| [1](../../validation/known-issues/category-1/index.md) | [Output](../../validation/known-issues/category-1/index.md) | Where XFOIL's printed and plotted output differs from the solver's state, or is too coarse to compare against. |
| [2](../../validation/known-issues/category-2/index.md) | [Persistent state](../../validation/known-issues/category-2/index.md) | Variables XFOIL leaves unset, stale or carried over between calls, which a translation has to reproduce. |
| [3](../../validation/known-issues/category-3/index.md) | [Drela's fixes](../../validation/known-issues/category-3/index.md) | Fixes and defects Drela marked in the shipped source, and how yFoil treats each one. |
| [4](../../validation/known-issues/category-4/index.md) | [Silent recovery](../../validation/known-issues/category-4/index.md) | Places where XFOIL corrects something or gives up silently, and carries on. |
| [5](../../validation/known-issues/category-5/index.md) | [Dead code](../../validation/known-issues/category-5/index.md) | Code in XFOIL 6.99 that the analysis path never runs. |
| [6](../../validation/known-issues/category-6/index.md) | [Geometry](../../validation/known-issues/category-6/index.md) | Where XFOIL's geometry handling differs from what its names, or the NACA definitions, suggest. |
| [7](../../validation/known-issues/category-7/index.md) | [Numerics](../../validation/known-issues/category-7/index.md) | Numerical details a reproduction has to copy exactly, and behaviours they produce. |
| [8](../../validation/known-issues/category-8/index.md) | [Open items](../../validation/known-issues/category-8/index.md) | Known gaps in yFoil's own tests, fixtures and studies. |
| [9](../../validation/known-issues/category-9/index.md) | [Pathological cases](../../validation/known-issues/category-9/index.md) | Runs whose reference is pathologically slow or hangs. |
