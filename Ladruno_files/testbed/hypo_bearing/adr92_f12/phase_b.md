> [!failure] **PHASE (b)/(c) — REFUTED, and not narrowly.**
> On the CP1/ADR-95 bearing deck scheme 2 committed **11 steps to `s/B = 4e-5` in
> 1347 s**, against the scheme-1 baseline's **51 steps to `s/B = 0.019` in 1267 s**
> — **475× shallower for the same wall clock**. Its step size collapsed to 25× the
> subdivision floor (`ds = 0.005 mm` against the baseline's `2.5 mm`, `free =
> False`), **every one of its 11 committed steps needed the relaxed rung 3**
> (100 % past rung 1, against CP1's 61–83 %), and it spent 4 of 80 subdivisions
> doing it. **The 1 % load–settlement bar could not be evaluated: there is no
> overlap in `s/B` between the arms at all.**

**Deck.** `sanisand_tau0_band.py`, leg `h1.0_e0.6944`, purpose-sized domain
`--xlim 10 --zbot 8` (624 DOF, 84 hexes — the `adr92_deck/x10z8` configuration),
`--maxsubsteps 20000`, `NormUnbalance @ 1e-5·γV`, `LoadControl` (stock, no
`-tangentPredictor`), 1200 s budget per arm, arms run **sequentially** so the wall
column means something. `INT_SCHEME` is a module constant with no CLI or env seam
(`:348`, read at call time at `:696` and `:1090`); `f12_bvp.py` loads the driver as
a module and sets it, so **no repo file was edited** and every one of the driver's
own controls, assertions and JSON provenance still ran. All three arms cleared
the geostatic controls identically (`resultant 4.44e-16`, `patch 1.89e-14`,
`η/M_c 0.6425`, `OutsideBounding 0`, `CLAMPING 0`).

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

**Reading.**

1. **The two scheme-1 arms are a working control.** `TanType 0` and `TanType 2`
   agree to **0.8–1.7 %** at every matched `s/B` and reproduce CP1's ladder rates
   (61.5 % / 70.6 % past rung 1, 85.4 % / 88.2 % failed-rung iteration share
   against CP1's 61–83 % / 89–93 %). So the reduced deck is representative and the
   instrument is sound. *(Aside: on this 84-element deck `TanType 0` is not the
   trap the driver's comment describes — it is marginally **faster** to every
   checkpoint. The driver's measurement was at `h0 = 0.25`, where the linear
   solve is not negligible; it is not contradicted, only bounded.)*
2. **Scheme 2 never leaves the subdivision floor.** `relaxed = 1` on every
   committed row from the first; `ds` fell from `0.02 mm` to `0.005 mm` and stayed
   there. This is §4's mechanism on a real BVP: the global Newton's trial iterates
   are off-path, the CPPM cannot return them, it grinds its recursive-halving
   ladder, the rung fails, the controller halves `ds`, and the smaller step does
   not help because the *iterate*, not the step, is what the CPPM chokes on.
   122.5 s per committed step against the baseline's 23–25 s.
3. **`base_foot_mismatch` is 7.9e-2 on the scheme-2 arm against 2–3e-4 on the
   others.** Do **not** read this as an equilibrium defect: the arm is at
   `q = 2.29 kPa`, a factor ~260 below the others, so the same absolute residual
   is a far larger relative one. It is reported because the driver reports it, not
   as a finding.
4. **`-implex` was OFF in all three arms.** Whether scheme 2 is the better
   *companion* — the question §9.2 leaves open — is **not** answered here and
   cannot be: the companion runs at `commitState` on an increment that is already
   given, which is §2's replay regime (where scheme 2 wins), not this one.
