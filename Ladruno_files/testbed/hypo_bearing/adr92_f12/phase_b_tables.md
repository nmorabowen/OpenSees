### 7.1 Leg outcome (x10z8, `h1.0_e0.6944`, 624 DOF, 84 hexes, 1200 s budget each)

| arm | scheme | TanType | steps | mode | s/B reached | q at end (kPa) | wall s | s per step | nfail | nsub | nrelax | CLAMPING | OutsideBounding |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| IntScheme 1, TanType 0 | 1 | 0 | 52 | WALL | 0.02155 | 668.78 | 1209 | 23.24 | 50 | 1 | 15 | 0 | 0 |
| IntScheme 1, TanType 2 | 1 | 2 | 51 | WALL | 0.01905 | 593.02 | 1267 | 24.84 | 61 | 1 | 22 | 0 | 0 |
| IntScheme 2, TanType 2 | 2 | 2 | 11 | WALL | 0.00004 | 2.29 | 1347 | 122.50 | 34 | 4 | 11 | 0 | 0 |

### 7.2 Committed load-settlement at matched s/B

| s/B | IntScheme 1, TanType 0 | IntScheme 1, TanType 2 | IntScheme 2, TanType 2 |
|---|---|---|---|
| 0.001 | 34.19 | 34.76 | -- |
| 0.002 | 64.53 | 65.44 | -- |
| 0.005 | 153.55 | 155.17 | -- |
| 0.010 | 303.38 | 305.86 | -- |
| 0.015 | 458.34 | 462.52 | -- |
| 0.020 | 618.30 | -- | -- |

**scheme 2 vs scheme 1 (both TanType 2) at matched s/B — the 1 % bar:**

| s/B | q s1_T2 | q s2_T2 | rel. diff | inside 1 %? |
|---|---|---|---|---|
| 0.001 | 34.76 | -- | -- | (not reached by both) |
| 0.002 | 65.44 | -- | -- | (not reached by both) |
| 0.005 | 155.17 | -- | -- | (not reached by both) |
| 0.010 | 305.86 | -- | -- | (not reached by both) |
| 0.015 | 462.52 | -- | -- | (not reached by both) |
| 0.020 | -- | -- | -- | (not reached by both) |

### 7.3 Wall clock to matched s/B

| s/B | IntScheme 1, TanType 0 | IntScheme 1, TanType 2 | IntScheme 2, TanType 2 |
|---|---|---|---|
| 0.001 | 27 s | 44 s | not reached |
| 0.002 | 68 s | 102 s | not reached |
| 0.005 | 258 s | 278 s | not reached |
| 0.010 | 549 s | 633 s | not reached |
| 0.015 | 842 s | 883 s | not reached |
| 0.020 | 1150 s | not reached | not reached |

### 7.4 Ladder decomposition (`adr92_bvp_gate.py`'s own arithmetic)

| arm | steps | rung1 | rung2 | rung3 | past rung 1 % | failed-rung iteration share % |
|---|---|---|---|---|---|---|
| IntScheme 1, TanType 0 | 52 | 20 | 17 | 15 | 61.5 | 85.4 |
| IntScheme 1, TanType 2 | 51 | 15 | 14 | 22 | 70.6 | 88.2 |
| IntScheme 2, TanType 2 | 11 | 0 | 0 | 11 | 100.0 | 95.7 |
