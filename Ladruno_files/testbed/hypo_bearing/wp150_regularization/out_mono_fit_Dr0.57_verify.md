## C++ verification of `fit_Dr0.57` (nb 2.070, A0 0.960, nd 2.388, h0 3.237); build 6bd6905b3ff8e5e3bdcc51bf01b05264e2a04177

| test | max \|Δq\|/q_max, C++ R1 off vs exact oracle | R1 on vs off |
|---|---|---|
| TX_p10 | 1.79e-03 | 9.06e-05 |
| TX_p50 | 1.58e-03 | 1.93e-04 |
| TX_p150 | 1.44e-03 | 2.35e-04 |
| TX_p500 | 1.25e-03 | 1.58e-05 |
| PS_p10 | 1.89e-03 | 7.65e-05 |
| PS_p50 | 1.44e-03 | 5.21e-05 |
| PS_p150 | 1.04e-03 | 1.56e-05 |
| PS_p500 | 8.79e-04 | 1.20e-05 |

worst: C++ vs oracle 1.89e-03; R1 on vs off 2.35e-04

| TXu monotonic, e 0.6944, p0 100 (R1 on) | p′_min (kPa) at ε_a | p′ at \|ε_a\| 1 / 5 % |
|---|---|---|
| TXu_comp_p100 | 38.2 at 0.0015 | 558 / 4551 |
| TXu_ext_p100 | 60.7 at 0.0013 | 570 / 4209 |

| CTXu (e 0.6944, p0 100) | N at 5 % DA, R1 off | R1 on | p′_min off / on |
|---|---|---|---|
| CTXu_csr0.20 | 0.5 | 0.5 | 1.0 / 1.0 |
| CTXu_csr0.15 | 1.0 | 1.0 | 1.0 / 1.0 |
