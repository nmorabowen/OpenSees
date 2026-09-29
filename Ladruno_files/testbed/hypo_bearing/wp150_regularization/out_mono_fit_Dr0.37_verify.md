## C++ verification of `fit_Dr0.37` (nb 1.134, A0 0.892, nd 1.178, h0 4.628); build 6bd6905b3ff8e5e3bdcc51bf01b05264e2a04177

| test | max \|Δq\|/q_max, C++ R1 off vs exact oracle | R1 on vs off |
|---|---|---|
| TX_p10 | 1.54e-03 | 6.52e-06 |
| TX_p50 | 1.47e-03 | 3.41e-04 |
| TX_p150 | 1.23e-03 | 2.26e-04 |
| TX_p500 | 1.27e-03 | 1.57e-05 |
| PS_p10 | 1.82e-03 | 8.45e-05 |
| PS_p50 | 1.19e-03 | 9.44e-06 |
| PS_p150 | 1.09e-03 | 1.10e-05 |
| PS_p500 | 9.72e-04 | 1.45e-05 |

worst: C++ vs oracle 1.82e-03; R1 on vs off 3.41e-04

| TXu monotonic, e 0.6944, p0 100 (R1 on) | p′_min (kPa) at ε_a | p′ at \|ε_a\| 1 / 5 % |
|---|---|---|
| TXu_comp_p100 | 28.6 at 0.0017 | 326 / 3584 |
| TXu_ext_p100 | 52.5 at 0.0016 | 318 / 2835 |

| CTXu (e 0.6944, p0 100) | N at 5 % DA, R1 off | R1 on | p′_min off / on |
|---|---|---|---|
| CTXu_csr0.20 | 0.5 | 0.5 | 0.9 / 0.9 |
| CTXu_csr0.15 | 1.0 | 1.0 | 0.9 / 0.9 |
