---
title: "ADR 92 / P2-9 — control-informed extrapolation factor: the oracle lane"
project: Ladruno
type: measurement
status: "ORACLE LANE COMPLETE — variant D implemented and measured; one refuted premise, one PASS with a 60x margin; C++ lane NOT started"
priority: high
owner: nmora
related:
  - "[[_adr92_p2_9_control_informed_f_plan]]"
  - "[[92_ladruno_sanisand_implex_adr]]"
  - "[[_adr92_p0_oracle_results]]"
  - "[[_adr92_p2_direction_oracle]]"
  - "[[_adr93_seat_replay]]"
  - "[[LadrunoSANISAND_implex_guide]]"
tags: [adr, sanisand, implex, p2, oracle, measurement]
updated: 2026-09-07
---

# ADR 92 / P2-9 — the control-informed factor `f*` in the numpy oracle

> [!abstract] **Headline.** The operator works and the seat prediction passes with room to
> spare (`0.4625 → 0.00082`, the plan asked for `≤ 0.05`), but the run turned up a
> **failure mode the plan does not name**: `f*` is only as good as the `Δε` it is frozen
> on. Frozen on a *bad* first iterate it collapses to `f* ≈ 0` and makes the path error
> **100x worse than today's `f`**; recomputed on the converged `Δε` it is `f* ≈ 0.99` and
> is **1.0–2.1x better than today's `f`** on the same rows. The `f* ≤ min(err(0),
> err(f_max))` result is an *identity at one trial*, not evidence about a path — the
> distinction is the whole finding.

## 1. How variant D is implemented (10 lines)

`Implex.FORMS` gains `"D"`. It keeps form A's *direction* — the committed `Δε_p(n)` — and
replaces A's *degree*:

1. `Implex.companion(deps)` runs the implicit return at the trial (`Sanisand.integrate`
   from the committed state) and restores the counters, so it is side-effect free. This
   is what `-implexControl` already computes; without control there is no companion.
2. `Implex.control_factor(deps, f_max)` forms `B = Ce:Δε_p(n)`, `BB = B:B`.
3. If `BB ≤ 0` it returns `f* = 0`, `A:B = 0` **without probing the companion** — the
   extrapolation term is identically zero for every `f`, so the factor is immaterial.
4. If the companion raises `Abandoned` it returns `f_max` (today's behaviour) and counts
   it in `d_probe_fail`. Falling back on `0` would silently turn IMPL-EX into an *elastic
   predictor* wherever the return map struggles — the opposite of the intent.
5. Otherwise `A = σ_n + Ce:Δε − σ_impl`, `f* = clamp((A:B)/(B:B), 0, f_max)`.
6. `A:B` is `dd_contr` — the **contravariant** double contraction, shear counted twice.
   Both operands are stress-like and `norm_contr(A − f B)² == dd_contr(A−fB, A−fB)`, so
   `f*` minimises **exactly** the quantity `implexError` measures. (Using `to_cov`/a
   strain-side product here would be wrong by a factor of 4 on every shear term.)
7. `extrapolate()` computes `f*` on the **first** `Δε` of the step and caches it in
   `self.f_star`; every later iterate of the same step reuses it, so the step stays
   linear in the frozen `Ce`. `commit()` clears it.
8. `Implex(..., freeze=False)` (driver kw `implex_freeze=False`) is the plan's priced
   alternative: recompute `f*` at every iterate. Diagnostic only — it costs linearity.
9. Census, taken once per step at `commit()` on the factor actually used: `d_backoff`
   (`f* < 0.5 f_max`), `d_zero`, `d_full`, `d_nullB`, `d_probe_fail`; `f_hist` records
   `(f_max, f*, A:B, B:B)` per step and `hist[-1]["f_used"]` is `implexDetail[5]`.
10. **Default OFF.** `Implex` still defaults to `form="A"`; `gate_GD` is opt-in and is
    *not* in `--gate all`, so the P0 memo's reproduction command is untouched.

Files: `Ladruno_implementation/adr92_p0_oracle/sanisand_implex_oracle.py`
(`Implex` docstring + `__init__`, `companion`, `control_factor`, the `"D"` branch of
`extrapolate`, the census in `commit`, `implex_freeze` on both drivers, `READ_GD4` +
`gate_GD` with `_gd_operator_controls` / `_gd_elastic_identity` / `_gd_reversal`), and
`Ladruno_implementation/adr93_p0/seat_replay.py` (`f_star`, `p29_seat`, `p29_ladder`,
`p29_forward`, `--only p29`).

## 2. Run commands

The G0/G2 probe CSVs are **not in the repo** (`adr92_p0_oracle/.gitignore` holds `data/`),
so they were regenerated first with the *installed* build — the same
`e95a1c74f7e15d7de8655eeeb004d7f34d81d512` the P0 memo cites. No C++ was built.

```bash
# 0. regenerate the P0 seeds (installed build e95a1c74f; ~3 s each)
P="C:/Users/nmb/venv/opensees_env/Scripts/python.exe"
$P Ladruno_implementation/adr92_p0_oracle/probe_binary_triaxial.py --p0 100            --nstep 40 --ez-max 0.02
$P Ladruno_implementation/adr92_p0_oracle/probe_binary_triaxial.py --p0 5              --nstep 40 --ez-max 0.02
$P Ladruno_implementation/adr92_p0_oracle/probe_binary_triaxial.py --p0 100 --e-init 0.6 --nstep 40 --ez-max 0.02

# 1. the untouched gates (byte-identity controls)
python3.12 -u Ladruno_implementation/adr92_p0_oracle/sanisand_implex_oracle.py --gate G0
python3.12 -u Ladruno_implementation/adr92_p0_oracle/sanisand_implex_oracle.py --gate G2

# 2. variant D
python3.12 -u Ladruno_implementation/adr92_p0_oracle/sanisand_implex_oracle.py --gate GD

# 3. the seat
python3.12 -u Ladruno_implementation/adr93_p0/seat_replay.py --only p29
```

## 3. G0 — the implicit gate, and the byte-identity controls

`gate_G0` PASSES on all three regenerated seeds, at the P0 memo's own precision:

| file | n | σ | e | α | z | εE | raw ≤ 1e-8 |
|---|---|---|---|---|---|---|---|
| `tx_p100_e0.6944_s1_n40.csv` | 40 | 4.19e-14 | 0 | 5.24e-15 | 1.10e-13 | 3.94e-14 | yes |
| `tx_p100_e0.6_s1_n40.csv` | 40 | 6.20e-14 | 0 | 9.06e-15 | 1.06e-14 | 4.42e-14 | yes |
| `tx_p5_e0.6944_s1_n40.csv` | 40 | 3.94e-13 | 0 | 3.55e-15 | 2.17e-14 | 1.75e-13 | yes |

`--gate G0` and `--gate G2` are **byte-identical before and after** the variant-D commit
(`diff` on the captured logs returns nothing; the G2 sweep is a full 21-row run at
`p0 = 100/5` and `e = 0.60`). That is the honest reading of the plan's "G0 row": G0
replays the binary's *implicit* path and contains no IMPL-EX at all, so variant D cannot
move it — the row is a **no-footprint control**, not a test of `f*`. The two limits the
row actually names are tested directly in GD.1/GD.2 below.

### GD.1 — the operator's limits, exactly

`σ_impl` is chosen so that `A = c·B` for a ladder of `c`; the operator must return
`clamp(c, 0, f_max)`. `|B| = 4.3364 kPa`, `B:B = 1.8804e+01`.

| c | f_max | f* returned | f* wanted | \|err\| |
|---|---|---|---|---|
| −2.00 | 1.0 | 0.00000000 | 0.00 | 0 |
| −0.25 | 1.0 | 0.00000000 | 0.00 | 0 |
| 0.00 | 1.0 | 0.00000000 | 0.00 | 0 |
| 0.25 | 1.0 | 0.25000000 | 0.25 | 1.6e-15 |
| 0.50 | 1.0 | 0.50000000 | 0.50 | 2.2e-16 |
| 1.00 | 1.0 | 1.00000000 | 1.00 | 0 |
| 2.00 | 1.0 | 1.00000000 | 1.00 | 0 |
| 0.75 | 0.5 | 0.50000000 | 0.50 | 0 |
| 0.30 | 0.5 | 0.30000000 | 0.30 | 5.0e-16 |

`B:B = 0` → `f* = 0`, `A:B = 0`, and the companion is probed **0 times** (no division,
no wasted return map). **PASS.**

### GD.2 — `B:B = 0` in a real run

An isotropic elastic unloading path (inside the yield cone): `max|σ_A − σ_D| = 0.0e+00`
kPa at n = 20 and n = 40, `B:B = 0` on 20/20 and 40/40 steps. **Bitwise identical.
PASS.**

### GD.3 — direction test (`drive_prescribed`, ONE trial per step)

| path | plastic steps | `f* = f_max` | `f* < 1e-3` | mean `f*` |
|---|---|---|---|---|
| monotone compression | 79 | 71 | 0 | 1.0000 |
| compression + reversal | 79 | 74 | 1 | 0.9873 |

The reversal neighbourhood — the turn is at step 41:

| step | f_max | f* | A:B | B:B |
|---|---|---|---|---|
| 39 | 1.000 | 0.999961 | +5.1588e+02 | 5.1590e+02 |
| 40 | 1.000 | 0.999949 | +5.2357e+02 | 5.2360e+02 |
| **41** | 1.000 | **0.000000** | **−1.7233e+02** | 5.3137e+02 |
| 42 | 1.000 | 1.000000 | +1.0356e+02 | 5.5610e+01 |
| 43 | 1.000 | 1.000000 | +2.2087e+02 | 1.8918e+02 |

Exactly the section-3 test the plan asks the C++ for: **`f* = 0` at the one step where the
history points the wrong way (`A:B < 0`), `f* = f_max` everywhere the history is right.**
This is the operator with no predictor pathology, and it is clean. **PASS.**

## 4. GD.4 — the G2 sweep at p0 = 100 and p0 = 5 (the surprise)

`TOT` = error vs the implicit reference at `dεz = 1e-5`; `ie mean` = mean `implexError`.
`D` freezes `f*` on the first iterate (the plan's operator); `Dr` recomputes it every
iterate, so the factor it commits was formed on the **converged** `Δε`.

**T1, p0 = 100 kPa**

| dεz | N | TOT A | **TOT D** | **TOT Dr** | TOT impl | ieA mean | ieD mean | ieDr mean | f* D | f* Dr | zeroD | zeroDr |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1.0e-4 | 200 | 1.511e-3 | **1.535e-1** | **1.572e-3** | 1.538e-3 | 5.186e-4 | 6.673e-2 | 5.150e-4 | 0.0353 | 0.9999 | 157 | 0 |
| 2.0e-4 | 100 | 6.936e-3 | 3.253e-1 | 6.930e-3 | 6.863e-3 | 2.314e-3 | 1.276e-1 | 2.053e-3 | 0.0210 | 0.9892 | 76 | 1 |
| 5.0e-4 | 40 | 5.123e-2 | 1.013e+0 | 5.120e-2 | 4.662e-2 | 1.540e-2 | 2.616e-1 | 1.095e-2 | 0.0115 | 0.9705 | 29 | 1 |
| 1.0e-3 | 20 | 2.319e-1 | 4.245e+0 | 2.905e-1 | 1.945e-1 | 6.453e-2 | 2.690e-1 | 3.132e-2 | 0.0000 | 0.9348 | 14 | 1 |
| 2.0e-3 | 10 | 4.056e+1 | 3.247e+1 | 2.752e+1 | 7.854e-1 | 3.151e-1 | 1.509e-1 | 1.615e-1 | 0.7156 | 0.9498 | 1 | 0 |
| 5.0e-3 | 4 | 1.589e+1 | 2.604e+1 | 2.530e+1 | 7.244e-1 | 3.572e-1 | 2.664e-2 | 2.688e-2 | 0.0000 | 0.3333 | 2 | 1 |
| 1.0e-2 | 2 | 1.775e+1 | 1.775e+1 | 1.775e+1 | 1.775e+1 | 4.56e-8 | 0 | 0 | 0.0000 | 0.0000 | 1 | 1 |

**T2, p0 = 5 kPa**

| dεz | N | TOT A | **TOT D** | **TOT Dr** | TOT impl | ieA mean | ieD mean | ieDr mean | f* D | f* Dr | zeroD | zeroDr |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1.0e-4 | 200 | 2.441e-3 | **1.103e+0** | **2.473e-3** | 2.213e-3 | 2.748e-3 | 1.831e-1 | 2.037e-3 | 0.0000 | 0.9942 | 176 | 1 |
| 2.0e-4 | 100 | 1.438e-2 | 2.920e+0 | 1.423e-2 | 1.051e-2 | 1.104e-2 | 2.546e-1 | 5.911e-3 | 0.0000 | 0.9876 | 86 | 1 |

> [!warning] The T2 rows at `dεz ≥ 5e-4` were **still running when this memo was written**
> and are not reported. That is the regime the P0 memo already calls "unusable from
> `5e-4` at `p0 = 5`" (hundreds of substeps per step), and the `Dr` arm re-probes the
> companion at every secant iterate on top of it. The two rows above are the campaign's
> own increments (`1e-4` = h 1.0 m nominal, `2e-4` = h 0.5 m nominal) and carry the
> verdict; re-run `--gate GD` with a longer budget to complete the block.

Same story at the corner, sharper: frozen-on-a-bad-predictor `D` is **450x** worse than A
(`1.103e+0` vs `2.441e-3`) with `f*` mean **0.0000** and `f* = 0` on 176/200 steps;
recomputed `Dr` is `2.473e-3` — 1.3 % off A on the path error and **1.3–1.9x better than A
on the control error** (2.037e-3 vs 2.748e-3; 5.911e-3 vs 1.104e-2) with `f*` mean 0.994.

**Why D collapses here, and why it is not the operator's fault.** `drive_triaxial` is
*mixed*-control: `solve_lateral` secant-iterates `dεxx` so that `σxx = p0`, and its FIRST
iterate is `dεxx = 0` — an **oedometric** increment, not the drained-triaxial one. Under a
pure axial increment at constant lateral strain `p` rises faster than `q`, `η` DROPS, and
the step is nearly elastic; `σ_impl` at that trial carries almost no plastic strain, so
`A ≈ 0` and `f* ≈ 0`. And `f = 0` in IMPL-EX is **not** the implicit answer — it is the
*elastic predictor*, which is far stiffer, so the lateral solve then converges on the
wrong strain and the path error explodes (`1.5e-1` against `1.5e-3`). Column `Dr` removes
exactly one thing — the bad predictor — and every symptom goes with it: `f*` returns to
0.93–1.00, `zeroD` 157 → `zeroDr` 0, and `TOT Dr` lands on `TOT A` to within 4 % while
`ieDr mean` is **1.0–2.1x better than A** on every resolved row (5.150e-4 vs 5.186e-4;
2.053e-3 vs 2.314e-3; 1.095e-2 vs 1.540e-2; 3.132e-2 vs 6.453e-2).

The last two rows of each block (`dεz ≥ 5e-3`) are past the driver's own breakdown —
`TOT impl` is `O(1)` there and the run prints `!! CONTAMINATED RUN` — so read them as
breakdown, not as data.

## 5. The seat (ADR-93 step 331 → 332): the registered prediction

State: element 4095 GP 8, `p = 56.70 kPa`, `|σ| = 115.75`, `f = 1.0` at step 332,
binary's committed `err = 0.4625`, `den = 121.683 kPa`. The section-4 identity
(`σ~ − σ_impl = Ce:(Δε_p(n+1) − f Δε_p(n))`, residual `1.8e-12 kPa`) makes `A` and `B`
exact, and the section-4 calibration supplies the seat's own `Δε_p(n)`.

```
|A| = |σ_n + Ce:Δε − σ_impl| =  3.4992 kPa      ( = |Ce:Δε_p(n+1)| )
|B| = |Ce:Δε_p(n)|           = 59.7758 kPa      ( 17.1x |A| — the stale history )
A:B = +2.0908e+02   B:B = 3.5732e+03   cos(A,B) = +0.9996
f_max = 1.0   ->   f* = clamp(A:B/B:B, 0, f_max) = 0.058515
```

| arm | f | err | x binary |
|---|---|---|---|
| (0) as shipped | 1.000000 | **0.4625** | 1.000 |
| **(v) P2-9 `f = f*`** | **0.058515** | **0.00082** | **0.002** |
| floor: no extrapolation | 0.000000 | 0.0288 | 0.062 |

**`cos(A, B) = 0.9996`** is the physical statement the seat has been missing: at the seat
the history's *direction* is right to within 1.6°; only its *magnitude* is wrong, by 17x.
That is precisely what a scalar factor can fix — and it is the independent reason
"trial-direction B" was rejected in `_adr92_p2_direction_oracle`.

### `f*` along the seat path (history magnitude swept)

Same state, same `Δε`, the history magnitude `t` swept; the seat sits at 20x.

| \|Δε(n)\| | /step332 | \|B\| | cos(A,B) | f* | err D | err f_max | err f=0 | substeps |
|---|---|---|---|---|---|---|---|---|
| 1.628e-5 | 0.25 | 0.238 | 0.9422 | 1.000000 | 0.0269 | 0.0269 | 0.0288 | 1 |
| 3.256e-5 | 0.50 | 0.784 | 0.9200 | 1.000000 | 0.0230 | 0.0230 | 0.0288 | 1 |
| 6.513e-5 | 1.00 | 3.499 | 1.0000 | 1.000000 | 0.0000 | 0.0000 | 0.0288 | 43 |
| 1.303e-4 | 2.00 | 8.809 | 0.9999 | 0.397192 | 0.0004 | 0.0436 | 0.0288 | 83 |
| 2.605e-4 | 4.00 | 20.676 | 0.9997 | 0.169194 | 0.0007 | 0.1412 | 0.0288 | 159 |
| 5.210e-4 | 8.00 | 46.616 | 0.9997 | 0.075038 | 0.0008 | 0.3544 | 0.0288 | 334 |
| 7.815e-4 | 12.00 | 73.183 | 0.9995 | 0.047793 | 0.0009 | 0.5727 | 0.0288 | 486 |
| 1.042e-3 | 16.00 | 100.096 | 0.9994 | 0.034939 | 0.0010 | 0.7939 | 0.0288 | 628 |
| 1.303e-3 | 20.00 | 127.169 | 0.9994 | 0.027499 | 0.0010 | 1.0163 | 0.0288 | 763 |
| 1.563e-3 | 24.00 | 154.321 | 0.9993 | 0.022659 | 0.0011 | 1.2395 | 0.0288 | 890 |

**`f* < 0.5 f_max` on 7/10 rows = 70 %** (the census slot the plan asks for). `f*` is
exactly `f_max` while the history matches the step (`≤ 1x`) and then tracks `1/t`;
`err D ≤ 0.0269` across four decades of history magnitude while `err f_max` runs to
`1.24`.

### A multi-step continuation from the seat (40 steps, ONE trial per step)

The census dumps only two committed rows at the seat, so the multi-step table is
generated: continue from state 331 along the step-332 strain direction at the step's own
magnitude.

| step | f* | err D | err f_max | err f=0 |
|---|---|---|---|---|
| 332 | 0.000000 | 2.875e-2 | 2.875e-2 | 2.875e-2 |
| 333 | 1.000000 | 1.507e-2 | 1.507e-2 | 4.328e-2 |
| 334 | 1.000000 | 4.175e-3 | 4.175e-3 | 4.674e-2 |
| 337 | 1.000000 | 8.913e-4 | 8.913e-4 | 4.908e-2 |
| 342 | 1.000000 | 2.104e-4 | 2.104e-4 | 4.819e-2 |
| 352 | 1.000000 | 5.904e-5 | 5.904e-5 | 4.442e-2 |
| 367 | 1.000000 | 1.676e-5 | 1.676e-5 | 3.914e-2 |
| 371 | 0.999984 | 2.528e-5 | 2.529e-5 | 3.791e-2 |

39 plastic steps, `f* < 0.5 f_max` on **0**, `f* = 0` on 1 (the first step, whose history
is empty), companion refusals 0. mean err D `1.3997e-3` vs f_max `1.4000e-3` vs f=0
`4.3652e-2`. **On a smooth path D is today's `f`** — it costs nothing where nothing is
wrong, which is the property the P2-2 zero-guard does not have.

## 6. Verdicts against the plan's section-2 table

**G0 rows — "byte-identical when `B·B = 0` and where `A ∥ B`". PASS, with a correction to
the row.** No G0 row moved: `--gate G0` and `--gate G2` are byte-identical before and
after, and G0's three seeds reproduce the binary at 4.2e-14 … 3.9e-13. But G0 is an
*implicit* replay and contains no IMPL-EX, so it cannot in principle test either limit —
the row as written is a no-footprint control. The limits themselves are PASS by direct
test: `B:B = 0` gives bitwise-identical stress on a real elastic path (GD.2, `max|σ_A −
σ_D| = 0.0e+00` on 60 steps) and never divides or probes the companion; `A = c·B` returns
`clamp(c, 0, f_max)` to ≤ 1.6e-15 over nine values of `c` including both clamps (GD.1).

**Seat replay step 331 — "error 0.46 → ≤ 0.05 (the `f = 0` value was 0.029); refuted if >
0.1". PASS, by 60x.** `f* = 0.058515` gives `err = 0.00082`: 560x below the binary's
0.4625, 35x below the `f = 0` floor of 0.0288, and 61x below the plan's 0.05 threshold.
**But state the caveat plainly:** at a single trial `f*` is by construction the minimiser
of `|A − f B|` over `[0, f_max]`, and both `0` and `f_max` lie in that interval, so
`err(f*) ≤ min(err(0), err(f_max))` is an **identity, not a measurement**. What was
measured is the *margin*: 0.00082 against the 0.0288 floor, i.e. `sin∠(A,B) = 0.0287` —
the seat's history is 17x too long but only 1.6° off axis.

**A new row the plan does not have — the predictor. PARTIAL / a live risk.** On the G2
drained-triaxial sweep, `f*` frozen on the first iterate is **100x worse than today's `f`**
(`TOT` 1.535e-1 vs 1.511e-3 at `dεz = 1e-4`, `f*` mean 0.0353, `f* = 0` on 157 of 200
steps), and recomputing it on the converged `Δε` recovers everything and then some (`TOT`
1.572e-3, `f*` mean 0.9999, and `implexError` 1.0–2.1x *better* than A on every resolved
row). The oracle's first iterate is pathological by construction (`dεxx = 0`, an
oedometric probe that unloads `η`), which a real global Newton's tangent predictor is not
— so this is not a refutation of P2-9. It is a **precondition on it**: freezing `f*` is
safe only if the first iterate's `Δε` is a *predictor* of the step, and the C++ must be
able to show that or must not freeze. See section 7.

**Fork R3 arm / Esmeralda dense + loose arms.** Out of the oracle lane — they need the
`-implexFactor control` binary. Not attempted; no C++ was built.

## 7. What the C++ implementer must carry over

1. **The inner product is `dd_contr`, contravariant.** `A:B = A11B11 + A22B22 + A33B33 +
   2(A12B12 + A23B23 + A13B13)`. Both operands are stress-like. In `ManzariDafalias`'s
   storage this is the same contraction `GetNorm_Contr` squares. A covariant/engineering
   product here is wrong by 4x on every shear term and would bias `f*` on any
   non-triaxial path — which is every Gauss point in a footing.
2. **Guard `B:B = 0` before the division AND before the companion.** With no committed
   plastic history the extrapolation term is identically zero for every `f`, so return
   `f* = 0` (or `f_max`; they are the same stress) and *skip the companion probe* — that
   is a whole return map saved on every elastic point.
3. **Decide the companion-refusal fallback explicitly, and make it `f_max`.** If the
   companion return map refuses at the trial, falling back on `f* = 0` turns IMPL-EX into
   an elastic predictor exactly where the material is hardest — the wrong direction. Fall
   back on today's `f` and count it (`d_probe_fail`); the control's own refusal path then
   handles the step.
4. **`σ_impl` must be the companion at the SAME `Δε` that `f*` is applied to.** The oracle
   makes this exact: `σ~ − σ_impl = Ce:(Δε_p(n+1) − f Δε_p(n))` holds to 1.8e-12 kPa
   because `ManzariDafalias` holds `G, K` at their committed values through the substep
   loop. If `Ce` were refreshed inside the step the identity — and `f*` — would drift.
5. **Freezing is the risk, not the operator.** Add the `implexDetail[5]` census
   (`f* < 0.5 f_max`) *and* a diagnostic that reports `f*` recomputed at the converged
   iterate versus the frozen one on a first campaign. If the two disagree by more than a
   few percent on a real deck, freeze on the *predictor's* `Δε` is not safe there and
   `-implexFactor control` should recompute (at the cost of step linearity) or the
   fallback is P2-8's fixed threshold, as the plan's own decision rule already allows.
6. **Sign convention.** Everything above is in the material's internal (compression-
   positive) convention, the one `getState` dumps; `LadrunoSANISAND3D` negates on the way
   out. `f*` is a scalar ratio and is convention-blind, but `A` and `B` must be formed in
   the same one.
7. **Cost.** One extra implicit return per *step* (not per iterate) when frozen, on top of
   the control's own companion — which under `-implexControl` is *already computed at the
   trial*, so if the code reuses it, `f*` is free apart from two contractions.

## 8. Reproduction footprint

- `Ladruno_implementation/adr92_p0_oracle/sanisand_implex_oracle.py` — variant D, opt-in
  `--gate GD`. `--gate all` byte-identical to the P0 memo's run.
- `Ladruno_implementation/adr93_p0/seat_replay.py` — `--only p29`, sections 6–8. Sections
  1–5 unchanged.
- `Ladruno_implementation/adr92_p0_oracle/data/` is gitignored; regenerate with the three
  `probe_binary_triaxial.py` commands in section 2 (installed build `e95a1c74f`).
- No `SRC/`, no `tests/`, no banner line, no C++ built.

## Log

- 2026-09-07 — oracle lane run. Variant D implemented; GD.1/GD.2/GD.3 PASS; seat
  prediction PASS at 0.00082; the freeze-on-a-bad-predictor collapse found and priced.
