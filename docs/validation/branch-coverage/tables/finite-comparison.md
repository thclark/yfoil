| case | calls | ref. iterations | outcome | parts at | XFOIL ill-conditioned at (its twins) | max \|ΔCL\| | max \|ΔCD\| | max \|ΔCM\| | max \|ΔRMSBL\| / \|ΔGAM\| | worst diff/floor | ref. 1-ULP spread |
|---|---:|---:|---|---|---|---:|---:|---:|---:|---:|---:|
| `kt_n60_inviscid_a4_cl05_cl8` | 3 | 24 | match | — | no call | 0 | — | 0 | 0 † | 0 | 2.8e-5 |
| `naca0012_n60_a2_re1e6_type3` | 1 | 9 | match | — | no call | 5.2e-15 | 1.4e-15 | 1.2e-15 | 1.5e-12 | 6.5e-1 | 1.4e-5 |
| `naca0012_n60_a8_re1e4_damp` | 1 | 15 | match | — | no call | 2.4e-14 | 2.8e-15 | 2.5e-16 | 2.4e-10 | 5.3e-1 | 5.1e-8 |
| `naca23012_n60_a4_re1e6_xtr022_0001` | 1 | 6 | match | — | no call | 1.1e-16 | 1.7e-18 | 1.9e-17 | 1.4e-14 | 1.7e-2 | 1.2e-10 |
| `naca4412_n160_polar_down16_re1e6_iter100` | 34 | 244 | divergent | call 31 it. 21 | calls 31, 32, 33, 34 | 4.9e-12 | 3.5e-13 | 2.8e-13 | 2.7e2 | 2.1e0 | 4.6e23 |
| `naca4412_n60_a18_re3e6_iter40` | 1 | 32 | match | — | call 1 | 1.2e-14 | 1.4e-14 | 1.7e-15 | 3.0e-9 | 8.1e-2 | 6.7e-6 |
| `naca4412_n60_a20_re1e6_iter50` | 1 | 50 | divergent | call 1 it. 0 | call 1 | — | — | — | 2.4e0 | 3.6e-1 | 5.1e3 |
| `naca4412_n60_cl1_clm05_re1e6` | 2 | 26 | match | — | no call | 1.1e-9 | 7.8e-12 | 5.9e-11 | 1.7e-8 | 1.0e0 | 3.1e-8 |
| `naca64a010_n60_inviscid_type2_m03_am23_a0_a19` | 3 | 0 | match | — | call 1 | 0 | — | 0 | 0 † | 0 | 5.8e3 |
| `naca64a010_n60_type2_a0_re1e6` | 1 | 16 | MISMATCH | call 1 it. 16 | call 1 | — | — | — | 4.1e-4 | 4.2e-3 | 4.6e7 |
