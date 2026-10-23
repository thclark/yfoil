# 2.9 The closure log's `s2` at laminar BLVAR calls is not a value BLVAR uses — Out of scope

The reference's closure log (`xfoil_subroutine_log.dat`, cut into `tests/fixtures/subroutines/` by
`closures = true`) records BLVAR's inputs as they stand in COM2. At some laminar calls (`ityp = 1`,
the stations near the leading edge) `S2` holds a subnormal (about 3e-314) that changed, with every
output unchanged, when the instrumentation's dump list grew from 32 to 1000 entries and again when
patches 20 and 21 were added — a value that depends on memory layout, so not one XFOIL computed for that station. BLVAR(1) does not read `S2`
(the laminar closures take no shear stress), so the outputs and the test are unaffected; the input
is regenerated with the fixtures and `--verify` holds it for a given build. Noted 2026-10-06.
