---
icon: lucide/repeat
---

# DO loops

A Fortran `DO` loop is compiled with a control test, and gcov counts that test's two edges against the
`DO` line: *iterate again*, or *leave the loop by its count*. Those counts cannot tell an empty loop from
one that ran and finished (see [Counting branches](counting-branches.md#what-the-statuses-mean)), so
this page states directly, for every loop, whether its body ran on the set.

--8<-- "validation/branch-coverage/tables/do-loops.md"
