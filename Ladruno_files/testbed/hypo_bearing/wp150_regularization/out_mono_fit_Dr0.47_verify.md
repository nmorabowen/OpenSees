## C++ verification of `fit_Dr0.47` (nb 1.615, A0 0.940, nd 1.771, h0 3.751); build 6bd6905b3ff8e5e3bdcc51bf01b05264e2a04177

| test | max \|Δq\|/q_max, C++ R1 off vs exact oracle | R1 on vs off |
|---|---|---|
| TX_p10 | 1.66e-03 | 5.61e-04 |
| TX_p50 | 1.34e-03 | 4.68e-04 |
| TX_p150 | 1.24e-03 | 3.99e-05 |
| TX_p500 | 1.10e-03 | 1.54e-05 |
| PS_p10 | 1.88e-03 | 9.36e-05 |
| PS_p50 | 1.19e-03 | 3.29e-04 |
| PS_p150 | 1.07e-03 | 1.06e-05 |
| PS_p500 | 9.03e-04 | 1.16e-05 |

worst: C++ vs oracle 1.88e-03; R1 on vs off 5.61e-04

| TXu monotonic, e 0.6944, p0 100 (R1 on) | p′_min (kPa) at ε_a | p′ at \|ε_a\| 1 / 5 % |
|---|---|---|
| TXu_comp_p100 | 32.1 at 0.0016 | 454 / 4374 |
| TXu_ext_p100 | 55.6 at 0.0015 | 452 / 3752 |

| CTXu (e 0.6944, p0 100) | N at 5 % DA, R1 off | R1 on | p′_min off / on |
|---|---|---|---|
| CTXu_csr0.20 | 0.5 | 0.5 | 1.0 / 1.0 |
| CTXu_csr0.15 | 1.0 | 1.0 | 0.8 / 0.8 |
