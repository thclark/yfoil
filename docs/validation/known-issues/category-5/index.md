---
icon: lucide/skull
---

# 5. Dead code

Code in XFOIL 6.99 that the analysis path never runs.

| | Issue | Description | Handling |
|---|---|---|---|
| [5.1](5-1-subroutines-never-called.md) | [Subroutines never called](5-1-subroutines-never-called.md) | UECALC, DSSET and DIT are never called in 6.99. | Not translated |
| [5.2](5-2-limage-is-permanently-false.md) | [LIMAGE is permanently false](5-2-limage-is-permanently-false.md) | The ground-image toggle is commented out, so the image-aerofoil code is dead. | Out of scope |
| [5.3](5-3-stale-duplicate-source-files.md) | [Stale duplicate source files](5-3-stale-duplicate-source-files.md) | Several source files are stale copies that are not built. | Out of scope |
| [5.4](5-4-matyp-3-is-unreachable.md) | [MATYP = 3 is unreachable](5-4-matyp-3-is-unreachable.md) | Nothing sets MATYP = 3, so its branches are dead. | Replicated |
| [5.5](5-5-oper-damp-is-reachable-and-undocumented.md) | [OPER `DAMP` is reachable and undocumented](5-5-oper-damp-is-reachable-and-undocumented.md) | OPER's DAMP command is missing from the menu but reachable. | Replicated |
| [5.6](5-6-structurally-dead-branches.md) | [Structurally dead branches](5-6-structurally-dead-branches.md) | Branches that the data on the analysis path can never reach. | Out of scope |
