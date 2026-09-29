## C++ verification of `fit_Dr0.47_nd3.5` (nb 1.652, A0 0.692, nd 3.500, h0 3.500); build 6bd6905b3ff8e5e3bdcc51bf01b05264e2a04177

| test | max \|Δq\|/q_max, C++ R1 off vs exact oracle | R1 on vs off |
|---|---|---|
| TX_p10 | 1.66e-03 | 5.57e-04 |
| TX_p50 | 1.34e-03 | 3.60e-04 |
| TX_p150 | 1.25e-03 | 4.63e-04 |
| TX_p500 | 1.30e-03 | 1.53e-05 |
| PS_p10 | 1.90e-03 | 9.81e-05 |
| PS_p50 | 1.49e-03 | 7.45e-05 |
| PS_p150 | 1.31e-03 | 7.65e-05 |
| PS_p500 | 1.18e-03 | 2.14e-04 |

worst: C++ vs oracle 1.90e-03; R1 on vs off 5.57e-04

| TXu monotonic, e 0.6944, p0 100 (R1 on) | p′_min (kPa) at ε_a | p′ at \|ε_a\| 1 / 5 % |
|---|---|---|
| TXu_comp_p100 | 63.1 at 0.0014 | 508 / 3880 |
| TXu_ext_p100 | 79.0 at 0.0011 | 513 / 3610 |

| CTXu (e 0.6944, p0 100) | N at 5 % DA, R1 off | R1 on | p′_min off / on |
|---|---|---|---|
| CTXu_csr0.20 | 1.0 | 1.0 | 1.0 / 1.0 |
| CTXu_csr0.15 | 1.5 | 1.5 | 0.9 / 0.9 |
