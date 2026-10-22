---
icon: lucide/git-pull-request
---

# Contributing

This tab is for people working on yFoil itself: how to find and fix a translation bug, how the
tests are organised, how the reports are regenerated, the project's conventions, and how to build
this site. It assumes the [Fidelity to XFOIL](../fidelity-to-xfoil/index.md) tab, which describes
what the tests establish; this one describes how to work within them.

## The pages of Contributing, in reading order

| Page | Read it when |
| --- | --- |
| [Finding a translation bug](finding-a-translation-bug.md) | a value differs from XFOIL's, or before changing any translated code |
| [The test suite](the-test-suite.md) | you are writing a test: the six purposes, the rules every test follows, the layout |
| [Regenerating the reports](regenerating-the-reports.md) | a generated page is stale, or you changed what it reads |
| [Conventions](conventions/index.md) | you are naming something, writing a figure, committing, opening a pull request or filing an issue |
| [Building the documentation](building-the-documentation.md) | you changed a page and want to see it |

## Two rules before anything else

yFoil is a **translation**, and two rules follow from that, without exception:

- **Never modify the algorithm.** When something does not work, the temptation is to "fix" it —
  a relaxation factor, an iteration limit, a clamp, a simplified formula, a coupling term disabled
  "temporarily". That always masks the real bug and creates new ones. If XFOIL does X, yFoil does
  X, including XFOIL's own bugs; the one exception is listed with its evidence under
  [deliberate divergences](../fidelity-to-xfoil/xfoils-known-issues/deliberate-divergences.md).
- **Debug forwards, never from the end.** A final CD or CL says nothing about where a bug is;
  [finding a translation bug](finding-a-translation-bug.md) is the method.

The full statement of the project's rules is `CLAUDE.md` at the repository root.
