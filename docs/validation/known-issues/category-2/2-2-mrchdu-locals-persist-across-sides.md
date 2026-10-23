# 2.2 MRCHDU locals persist across sides — Replicated

`AMI, SENS/SENNEW, UEREF/HKREF, CTE/TTE/DTE` are Fortran locals never re-initialised between side 1
and side 2, and `COM1/COM2/XT` are COMMON (`src/bl/mrchdu.rs:64-65`). A per-side reset diverges from
XFOIL. The same persistence in MRCHUE is trace-faithful only: side 2's similarity station overwrites
COM1 before anything reads it (`src/bl/mrchue.rs:54-56`).
