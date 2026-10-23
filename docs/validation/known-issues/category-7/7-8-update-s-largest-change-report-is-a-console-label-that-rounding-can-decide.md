# 7.8 UPDATE's largest-change report is a console label that rounding can decide — Replicated

UPDATE records the single largest change it is about to make — its size RMXBL, the variable
VMXBL (`n`/`C`, `T`, `D`, `U`) and the station IMXBL, ISMXBL (`xbl.f:1433-1471`) — and VISCAL only
prints them (`xoper.f:3001`, `:3003`); nothing reads them back. Convergence is `RMSBL < EPS1`, the
root-mean-square change over every variable at every station (`xoper.f:3009`), and the relaxation
factor comes from the changes themselves. Two properties make the label a poor quantity to
compare: at the similarity station the relative changes of θ and δ* are mathematically equal, so
which one is reported is decided by the last bit of rounding; and for Ue the stored size is the raw
change `DUEDG` (`:1468`) while the candidates are ranked by the normalised `DN4`, so the reported
size changes units with the winner. **yFoil** replicates both (`update.rs`, `residual_max`). The
step tests compare the size's magnitude and the solution arrays, not the label; the label is
compared exactly on the well-conditioned reference case (`tests/execution/update.rs`). Not
corrected: yFoil outputs every distribution in full, so nothing depends on the console label.
