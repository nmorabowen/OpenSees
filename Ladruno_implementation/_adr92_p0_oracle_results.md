---
title: "ADR 92 / P0 — the single-Gauss-point IMPL-EX oracle: results"
project: Ladruno
type: results
status: "P0 COMPLETE — G0 PASS, G1 PASS, G2 measured, G3 measured (finding), G4 decides D1 = A, G5 PASS; D3 REVERSED"
priority: high
owner: nmora
related:
  - "[[92_ladruno_sanisand_implex_adr]]"
  - "[[_adr92_sanisand_implex_scoping]]"
  - "[[_adr90_tau0_qu_band]]"
tags: [adr, implex, sanisand, oracle, measurement]
updated: 2026-09-05
---

# ADR 92 / P0 — the oracle's answer

> [!success] **VERDICT**
> **D1 = A (extrapolate the plastic strain).** On the deck-default integrator the textbook
> `dGamma` form is **20x worse** at the calibrated confinement and **18x worse** at
> `p0 = 5 kPa`; under scheme 2 the two forms coincide to `2e-10 kPa`, which is the
> `mDGamma` finding of the scoping memo measured rather than argued.
>
> **D3 is REVERSED.** `BackwardEuler_CPPM` is not an implicit return where the campaign's
> problem lives: on the corner path its low-p branch — Newton disabled by a literal
> `errFlag = 0` (`:2264`) — sent **58–74 % of all calls to `explicit_integrator`**, i.e. to
> `ModifiedEuler`. The IMPL-EX companion defaults to **scheme 1 with `-maxSubsteps`**
> (#792 T1); scheme 2 is permitted and buys nothing at low `p`.
>
> **A finding the ADR did not carry: at the `p_min` floor the extrapolated stress goes
> through the floor into tension** — `min p = -1.37 kPa` at 40 steps, `-0.16` at 80, `-0.09`
> at 160, against a committed `p = +0.0101`. First order in the step, but O(1)–O(10) relative.
> **P1 must clamp `sigma~` with the code's own device** (`dev + p_min*I1`) or the element at
> the free surface receives tensile mean stress every iteration.
>
> **What P0 cannot say:** the seizure. `ModifiedEuler`'s force-accept at `dT_min` fired
> **zero** times on every prescribed path here (≤ 12 002 substeps / 200 steps). The
> `dT -> 1e-6` collapse GATE U measured needs the boundary-value problem's strain path at the
> corner; the "~100x" cost claim of ADR §2 stays **unmeasured until P3 / #792 T8**.

**Build:** every binary number is from `e95a1c74f7e15d7de8655eeeb004d7f34d81d512` read out of
`C:\Program Files\Ladruno\OpenSees\bin\opensees.pyd` by `ops.ladrunoBuild()`; interpreter
`C:\Users\nmb\venv\opensees_env\Scripts\python.exe` (3.12.10). The oracle is
`adr92_p0_oracle/sanisand_implex_oracle.py`, numpy 2.4.6 under `python3.12`; every table
below is `--gate <G>` of that script at the commit this memo lands in.

---

## 1. What was run

| piece | what |
|---|---|
| `probe_binary_triaxial.py` | the act's `probe_sanisand_tx.tcl` on `LadrunoSANISAND`, the D-L cell's 18 constants, `-Pmin 0.0101 -Presidual 0.0 -honorTolR 0`, `TanType 2`; dumps `eps(6) sig(6) state(26)` per step. Consolidation held at `LoadControl 0.0` — the Tcl original ramps 5 more steps and its "`p0 = 100`" is actually 150 kPa. |
| `sanisand_implex_oracle.py` | numpy transcription of `ManzariDafalias` at one Gauss point, **the code and not the paper**: `G` on `m_e_init`; the `D_factor` sigmoid; `ForwardEuler`'s zero-`r` **verbatim, scheme 5 only**; `TolE = 1e-4`; the force-accept at `dT_min`; `Stress_Correction` ON; `BackwardEuler_CPPM` with its low-p Newton disabled. Four integrators: `reference`, `implicit`, `implex_A`, `implex_B`. Compression positive inside; CSV columns negated on read. |
| G0 | oracle `implicit` replays the binary's own recorded strain sequence |
| G1 | error under increment halving, `N = 25…400` to `eps_z = 0.004`, reference at `6.25e-7` |
| G2 | error over the campaign's increment range, `1e-4 … 2e-3`, to `eps_z = 0.02` |
| G3 | a prescribed path onto the `p_min` floor with shear kept on, seeded from the `p0 = 5` state |
| G4 | A vs B, three states, plus a compression–reversal path |
| G5 | `d(sigma~)/d(eps)` against `Ce(p_n)` by central differences |
| instr | how often the force-accept and the scheme-2 fallback fire, per increment |
| controls | total-vs-incremental form; the `r` defect; A ≡ B under scheme 2 |

## 2. G0 — the oracle is not measuring itself

| binary path | n | `sigma` rel | `e` | `alpha` | `z` | `epsE` | verdict |
|---|---|---|---|---|---|---|---|
| `p100 e0.6944 s1` | 40 | **4.2e-14** | 0 | 5.2e-15 | 1.1e-13 | 3.9e-14 | exact |
| `p100 e0.6944 s1` | 400 | 2.9e-7 | 0 | 1.3e-8 | 4.0e-8 | 5.4e-6 | see control A |
| `p100 e0.60 s1` | 40 | **6.2e-14** | 0 | 9.1e-15 | 1.1e-14 | 4.4e-14 | exact |
| `p5 e0.6944 s1` | 40 | **3.9e-13** | 0 | 3.6e-15 | 2.2e-14 | 1.8e-13 | exact |
| `p1 e0.6944 s1` | 400 | **4.9e-13** | 0 | 2.0e-14 | 2.1e-14 | 1.5e-12 | exact |
| `p100 e0.6944 s2` | 40 | 9.6e-9 | 0 | 1.4e-8 | 1.3e-6 | 3.4e-7 | see control B |
| `p100 e0.6944 s2` | 400 | 1.0e-9 | 0 | 9.8e-8 | 9.1e-8 | 2.1e-8 | see control B |

**Control A** (the 400-step scheme-1 row): a `1e-15` relative perturbation of the seed moves
the oracle **from itself** by `2.6e-5`; the oracle-vs-binary gap is `5.4e-6`. The
substepper's accept/reject is discontinuous in the state, so a long path amplifies round-off
by ~10 orders; the gap is inside that. **Control B** (scheme 2): the code's own 19-component
residual evaluated at the binary's answer is `3.0e-8 / 4.2e-9`, at the oracle's `9.8e-9 /
1.5e-9` — both inside `NewtonIter2`'s stopping tolerance, the oracle's at least as small.
**G0 PASS.** Four of five deck-default paths are reproduced to round-off; nothing was tuned.

Two errors in `SOURCE_EXTRACTION.md` were found by G0 refusing to close and are corrected
there: `Stress_Correction` is **on** by default (`mStressCorrectionInUse = true`,
`:217/303/368/429`), and the scheme-2 low-p branch **always** integrates explicitly (`:2264`).

## 3. G1 — order of convergence

`EXTRAP` = extrapolated stress vs its own companion's answer (the IMPL-EX-specific error);
`TOTAL` = the run vs the fine reference (includes the substepper's own error).

**T1, `p0 = 100 kPa`, to `eps_z = 0.004`:**

| `d eps_z` | EXTRAP A | EXTRAP B (s1) | EXTRAP B (s2) | TOTAL A | TOTAL implicit | `implexError` mean |
|---|---|---|---|---|---|---|
| 1.6e-4 | 1.89e-4 | 2.20e-1 | 3.28e-2 | 1.90e-2 | 1.89e-2 | 7.0e-3 |
| 8.0e-5 | 5.03e-5 | 1.05e-1 | 1.20e-2 | 4.39e-3 | 4.34e-3 | 1.6e-3 |
| 4.0e-5 | 2.43e-5 | 5.20e-2 | 4.05e-3 | 2.93e-5 | 3.39e-5 | 3.2e-4 |
| 2.0e-5 | 5.30e-6 | 2.51e-2 | 2.08e-3 | 1.74e-5 | 1.57e-5 | 9.4e-5 |
| 1.0e-5 | 1.67e-6 | 1.25e-2 | 1.07e-3 | 2.52e-5 | 2.42e-5 | 2.5e-5 |
| **slope** | **1.69** (fine-3: 1.93) | 1.03 | 1.24 | — | — | **2.03** |

**T2, `p0 = 5 kPa`:** EXTRAP A `2.28e-2, 1.19e-3, 8.49e-5, 3.38e-5, 1.01e-4` (slope 2.08 —
non-monotone on the last point, which is the substepper's round-off floor showing through);
EXTRAP B (s1) `1.64 → 6.6e-2`, slope 1.15.

**Reading.** IMPL-EX-A's own error is first order or better (measured 1.7–2.1) and **`TOTAL A`
tracks `TOTAL implicit` to within 10 % at every increment** — at the campaign's step the
extrapolation adds essentially nothing to the error the substepper already carries. B on the
deck-default integrator is 10–100x worse and exactly first order; B on scheme 2 recovers
most of that, which is the `mDGamma` story. **G1 PASS.**

## 4. G2 — the campaign's increment

**T1, `p0 = 100 kPa`, to `eps_z = 0.02`, reference `1e-5`, terminal `q/p = 1.5500`:**

| `d eps_z` | N | EXTRAP A | EXTRAP B | `de/e` A | `d eta/eta` A | `implexError` A mean | `q/p` A | `q/p` impl | note |
|---|---|---|---|---|---|---|---|---|---|
| **1.0e-4** | 200 | **5.1e-5** | 1.5e-1 | 9.9e-6 | **9.2e-4** | 5.2e-4 | 1.5514 | 1.5516 | `h = 1.0 m` nominal |
| **2.0e-4** | 100 | **7.3e-5** | 3.1e-1 | 4.0e-5 | **4.7e-3** | 2.3e-3 | 1.5572 | 1.5571 | `h = 0.5 m` nominal |
| 5.0e-4 | 40 | 4.4e-3 | 9.2e-1 | 1.8e-4 | 3.3e-2 | 1.5e-2 | 1.6013 | 1.5968 | |
| 1.0e-3 | 20 | 3.2e-2 | 3.4 | 2.2e-4 | 1.3e-1 | 6.5e-2 | 1.7562 | 1.7270 | implicit itself 19 % off |
| 2.0e-3 | 10 | 2.3e+1 | 2.9e+1 | 7.3e-2 | 8.9e-1 | 3.2e-1 | 2.92 | 2.07 | **breakdown — the binary stalls here too** |

The last row is not an IMPL-EX failure: at `d eps_z = 2e-3` the implicit run is 79 % off its
own reference and the **binary's** single-element push at the same increment fails
`NormDispIncr` at step 2 (`probe … --nstep 100 --ez-max 0.20`). At `p0 = 1 kPa` the binary
stalls at `5e-4` and completes at `5e-5` — the campaign's low-confinement wall, at one Gauss
point. *(T2, `p0 = 5`, rows: `data/g2_rerun.log`, appended below when the rerun lands.)*

**Reading.** At the **nominal** campaign increment the extrapolation costs `5e-5` in stress
and `1e-3` in mobilised stress ratio — invisible against the substepper's own `1.5e-3`. The
README's caveat stands: the corner sees 10–100x the nominal strain increment, and at `1e-3`
the cost is 3 % in stress and 13 % in `eta`. `-implexControl` at `0.05` would refuse the
`5e-4` row and above.

## 5. G3 — the corner (the finding)

Path: seeded from the committed `p0 = 5 kPa` state; volumetric extension `-3e-3` reached by
`t = 0.5` and held, with deviatoric shear `6e-3·t` on throughout, so `p` collapses onto the
floor while the point keeps flowing. Reference at 10 240 steps.

| n | terminal `p` (ref / impl / A / B / s2) | **min `p` along path A / B** | A `‖dσ‖/‖σ‖` | ME force-accept | low-p clamps | s2 calls → low-p → **EXPLICIT** |
|---|---|---|---|---|---|---|
| 40 | 0.0101 / 0.0101 / 0.0101 / 0.0101 / 0.0101 | **−1.369 / −1.369** | 17.1 | 0 / 7289 | 19 | 86 → 58 → **58** |
| 80 | all 0.0101 | **−0.156 / −0.482** | 8.57 | 0 / 6159 | 57 | 196 → 121 → **121** |
| 160 | all 0.0101 | **−0.089 / −0.087** | 4.28 | 0 / 7201 | 102 | 256 → 189 → **189** |

Three things, in order of consequence:

1. **The extrapolated stress crosses the floor.** The committed state sits at `+0.0101 kPa`;
   nothing clamps `sigma~`, and `f·d_eps_p(n)` carries it to `−1.37 kPa` at the coarsest step.
   Halves with the step (first order) but is O(1)–O(10) relative at every step tried. The C++
   must apply the code's own floor device to `sigma~` — `dev(sigma~) + p_min*I1` — before it
   is handed to the element. **New P1 item; not in the request, not in the ADR draft.**
2. **Scheme 2 is `ModifiedEuler` at the corner.** 67 %, 62 %, 74 % of the companion's calls
   took the low-p branch, and every one of them integrated by `explicit_integrator`. The
   "only implicit return" argument for D3 does not survive contact with the place D3 was for.
3. **The force-accept never fired.** 0 of ~7 000 substeps on each path. The seizure GATE U
   measured is not a property of the material at a prescribed strain path; it is a property
   of the corner's strain path in the BVP. P0 cannot price it.

## 6. G4 — D1, head to head (`d eps_z = 1e-4`, 40 steps)

| state | A (s1) `‖dσ‖/‖σ‖` | A (s2) | B (s1) | B (s2) | A/B on s1 |
|---|---|---|---|---|---|
| `p0 = 100`, e 0.6944 | **6.7e-3** | 9.6e-3 | 1.4e-1 | 9.6e-3 | **21x** |
| `p0 = 5`, e 0.6944 | **5.0e-2** | 1.8e-2 | 8.9e-1 | 1.8e-2 | **18x** |
| `p0 = 100`, e 0.60 | **7.0e-3** | 1.3e-2 | 1.5e-1 | 1.3e-2 | **22x** |
| reversal `±0.004` | A 7.10e-2 (implicit itself 7.10e-2) | | B 9.6e-2 | | |

Under scheme 2, A and B are the same operator (`2.2e-10 kPa` apart — control). Under the
deck default they differ by **37 kPa on a 277 kPa stress**, because `mDGamma` is the last
substep's multiplier. On the reversal path A's error is the implicit's own to four figures:
the one-step lag of `alpha_in` costs nothing measurable here. **D1 = A.**

## 7. G5 — tangent identity

`max|d(sigma~)/d(eps) − Ce(p_n)| / max|Ce|` = **3.5e-11** (A), **5.7e-11** (B), central
differences at `1e-9`. **PASS.** The C++ keeps this only if `G`, `K` are read off the
**committed** stress; the base reads them there already (`:1008`, `:2223`).

## 8. Controls

- **Total-strain form** (the ADR's first draft): at a **zero** increment it returns
  `p = 99.50` against a committed `100.00` (0.5 %) and **`1.11` against `5.00` (78 %)**. The
  incremental form is not a refinement; it is the only correct one for a hypoelastic law.
- **`ForwardEuler`'s zero-`r`** (scheme 5 only): terminal `q/p` 0.946 (scheme 1, `r` correct)
  vs **1.188** (scheme 5). A 26 % error in mobilised strength, silent. Vanilla defect, owed a
  quirks row and a two-line fix; **does not reach the campaign.**
- **Instrumented counts on triaxial ramps**, `1e-4 … 2e-2`, three confinements: force-accept
  **0** everywhere; scheme-2 Newton failures and recursive halvings appear from `5e-4` up
  (2–10 per run) but **0 explicit fallbacks** — the fallback fires only through the low-p
  branch (G3), never through Newton failure on these paths.

## 9. What this licenses ADR 92 to say — and not

**MAY say:**
- D1 = plastic-strain extrapolation, by 18–22x on the deck-default integrator (§6).
- The IMPL-EX-specific error is first order or better (1.7–2.1 measured) and, at the nominal
  campaign increment, `5e-5` in stress / `1e-3` in `eta` — below the substepper's own (§3, §4).
- Scheme 2's low-p branch is explicit, so D3's rationale is void; the companion is scheme 1
  under `-maxSubsteps` (§5).
- The extrapolated stress needs the `p_min` clamp; without it the free-surface ring receives
  tensile mean stress (§5).
- The total-strain form is wrong by 78 % at `p0 = 5` (§8) — the negative control for P1.

**MUST NOT say:**
- Anything about the seizure, the "~100x", or wall-clock. Not measured; the force-accept
  never fired here.
- That the corner error is "small". It is first order and O(1)–O(10).
- That IMPL-EX reaches a plateau, a capacity, or a wall. Single Gauss point, prescribed
  strain, no BVP feedback.

## 10. The disclosure the act must print beside any IMPL-EX verdict

> This curve was produced with `LadrunoSANISAND -implex`. The stress handed to each element
> is extrapolated from the last committed step, not the converged state. At the nominal
> increment (`d eps ≈ 1e-4`) the extrapolation error is `5e-5` in stress and `1e-3` in
> mobilised stress ratio at `p0 = 100 kPa`; it grows to `4e-3 / 3e-2` at `5e-4` and
> `3e-2 / 0.13` at `1e-3`, and halves or better when the increment halves. At Gauss points on
> the `p_min` floor — the ring beside the footing — the extrapolated stress can leave the floor
> by O(1) kPa; that ring is the least trustworthy part of this field. Every equilibrium here
> is an equilibrium of the extrapolated stress; any limit point read off it must be confirmed
> on the implicit material to the last settlement the implicit solver reaches. The per-step
> `implexError` is printed beside this curve.

## 11. Reproducing

```bash
C:/Users/nmb/venv/opensees_env/Scripts/python.exe Ladruno_implementation/adr92_p0_oracle/probe_binary_triaxial.py --p0 100 --nstep 40 --ez-max 0.02
python3.12 Ladruno_implementation/adr92_p0_oracle/sanisand_implex_oracle.py --gate all
```

## Log

- 2026-09-05 — The P0 builder (Opus) ran G0 to a pass and was terminated by its session
  limit before the memo; Fable picked up the oracle as left, added the `exp` overflow guard,
  ran G1–G5 / instr / controls, and wrote this. Every number above is from this session's
  runs, not the builder's transcript.
