---
icon: lucide/circle-x
---

# Branches not taken

## Not taken by any measured run

No measured run takes these branches, and no code-reading argument says they cannot be taken. This is
a statement about the search for cases, not about the code. The reach note is the input the code
reading says would reach each one.

--8<-- "validation/branch-coverage/tables/branches-not-taken.md"

## Taken only in pathological cases

Pathological cases are copies of runs that hang in XFOIL's plot label or take hours
([known issues §9](../known-issues/category-9/index.md)), so they cannot
be tests. They are measured only on request.

--8<-- "validation/branch-coverage/tables/branches-pathological.md"
