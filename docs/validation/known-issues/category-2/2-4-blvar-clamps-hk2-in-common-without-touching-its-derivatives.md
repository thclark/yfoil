# 2.4 BLVAR clamps HK2 in COMMON without touching its derivatives — Replicated

`xblsys.f:798-799`: `HK2 = MAX(HK2, 1.00005)` in the wake, `MAX(HK2, 1.05)` otherwise, while
`HK2_U/T/D/MS` stay as BLKIN set them. The clamped value then feeds BLMID, the next station's HK1
(after `COM1 = COM2`) and MRCHUE's HTARG (`src/bl/system.rs:1171-1176`). The wake `US2` is likewise
clamped (`xblsys.f:833, 843`). The order of the repeated BLVAR calls in MRCHDU matters only because of
this in-place clamp (`src/bl/mrchdu.rs:385-387`).
