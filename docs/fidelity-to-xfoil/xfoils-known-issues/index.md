# XFOIL's known issues

XFOIL is not a clean specification. Its source has quirks, dated fixes, silent self-corrections,
uninitialised variables, dead branches, and behaviours that only show on particular runs — and a
faithful translation reproduces all of them. This section is where every one found during the
translation is recorded once, with how yFoil handles it: replicated, overcome, or out of scope.
Other pages that meet an issue link here rather than describe it again.

| Page | What it holds |
| --- | --- |
| [The issue register](the-issue-register.md) | Every known issue, by category, numbered (§1.1, §7.7 …), each with XFOIL's source, yFoil's handling and the test where one exists |
| [Deliberate divergences](deliberate-divergences.md) | The only place yFoil knowingly departs from XFOIL, and the standard of evidence for adding another |
| [Runs that show the issues](runs-that-show-the-issues.md) | Reference sweeps in which an issue arises — the panel-count dependence, the plot-label hang, the Kármán–Tsien pole — with yFoil compared point by point |

For a user the issues that matter most are the ones that make a returned number untrustworthy;
[whether to trust a result](../../using-yfoil/whether-to-trust-a-result.md) says how yFoil flags
those, and links back here for the detail. Each section of the register that has a test has one
module in `tests/known_issues/`.
