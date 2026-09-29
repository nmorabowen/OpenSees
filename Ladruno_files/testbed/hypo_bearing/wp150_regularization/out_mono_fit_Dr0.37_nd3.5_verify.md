## C++ verification of `fit_Dr0.37_nd3.5` (nb 1.174, A0 0.509, nd 3.500, h0 4.267); build 6bd6905b3ff8e5e3bdcc51bf01b05264e2a04177

| test | max \|Δq\|/q_max, C++ R1 off vs exact oracle | R1 on vs off |
|---|---|---|
| TX_p10 | 1.35e-03 | 9.42e-04 |
| TX_p50 | 1.28e-03 | 3.46e-04 |
| TX_p150 | 1.23e-03 | 1.52e-04 |
| TX_p500 | 1.14e-03 | 2.26e-04 |
| PS_p10 | 1.87e-03 | 1.02e-04 |
| PS_p50 | 1.46e-03 | 4.79e-05 |
| PS_p150 | 1.35e-03 | 4.64e-05 |
| PS_p500 | 1.24e-03 | 2.15e-05 |

worst: C++ vs oracle 1.87e-03; R1 on vs off 9.42e-04

| TXu monotonic, e 0.6944, p0 100 (R1 on) | p′_min (kPa) at ε_a | p′ at \|ε_a\| 1 / 5 % |
|---|---|---|
| TXu_comp_p100 | 73.7 at 0.0014 | 403 / 3316 |
| TXu_ext_p100 | 85.6 at 0.0011 | 394 / 2865 |

| CTXu (e 0.6944, p0 100) | N at 5 % DA, R1 off | R1 on | p′_min off / on |
|---|---|---|---|
| CTXu_csr0.20 | 1.5 | 1.5 | 0.9 / 0.9 |
| CTXu_csr0.15 | 2.0 | 2.0 | 0.8 / 0.8 |
