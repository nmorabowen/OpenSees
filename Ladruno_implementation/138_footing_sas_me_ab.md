---
title: "WP-138 — strip-footing A/B: ModifiedEuler (IntScheme 1) vs SAS-ME (IntScheme 129)"
project: Ladruno
type: measurement report
status: IN PROGRESS — local legs closed (control; A ModifiedEuler and B SAS-ME beb6d8333 both to s/B 0.01737); continuation on Esmeralda (E_A/E_B/E_C2/E_B16), results relayed by the orchestrator
related:
  - "[[_tims_2d_model_requests_2026-09-25]]"
  - "[[LadrunoSANISAND_implex_guide]]"
  - "[[LEDGER_quirks]]"
updated: 2026-09-28
---

# WP-138 — strip-footing A/B, ModifiedEuler vs SAS-ME (TIMs F18 end to end)

Pure Python decks and runs; no C++ change, no build. Deck, driver and every run
output: `Ladruno_files/testbed/footing_sas_me_ab/`. Nothing in the TIMs Workbench
was run or edited; the deck is built from the intake's §1 spec
(`_tims_2d_model_requests_2026-09-25.md`, on `origin/wp/127-sanisand-replay-counters`).

## 1. The deck

| item | value | source |
|---|---|---|
| geometry | plane strain, B = 1.5 m strip, domain 15B × 12B = 22.5 × 18 m, **full width** | intake §1 (the act saw asymmetric bands) |
| element | `LadrunoQuad -formulation bbar -type PlaneStrain`, 2 430 elements, 9 720 GP, 4 860 free DOF | intake §1 |
| mesh | 90 × 27 tensor mesh: 64 columns of B/8 over \|x\| ≤ 4B, 13 geometric columns per side to ±11.25 m; 12 rows of B/8 in the top 1.5B, 15 geometric rows to −18 m (first graded cell = B/8). Lower block numbered first (tags 1–1350, column-major), top block 1351–2430 | **reconstructed** (gap G1) |
| loads | buoyant self-weight γ' = 9.81 kN/m³ (ramped by `selfWeight`); 7.65 kPa on the surface outside the footprint; 18.4 kN/m on the footing's reference node | intake §1 / F10's transcription (G2) |
| K0 | ν* = 0.312885 → K0 = 0.4553 (Jaky, 33°), held for the whole run | attachments README |
| footing | `LadrunoKinematicCoupling` from a 3-DOF reference node at (0,0) to the 9 footprint nodes, `-dof 1 2` (rough), default penalty 1e12; reference u_x and θ fixed, u_y pushed by `sp` | G3 |
| BCs | base fixed; sides u_x = 0 | standard |
| solver | `system Pardiso`, `NormUnbalance` 1e-5 × the total applied vertical load (4 152.1 kN/m → 0.0415 kN); Newton (25 it) → NewtonLineSearch (40) → KrylovNewton (60, tol × 10) | intake §1, F10's ladder |
| controller | ds0 = 2e-5 m, ×2 after 6 good steps up to 1e-3 m, ÷2 on a failed ladder; stop when ds < 2e-7 m (FLOOR, "ds below floor") or on the wall-clock cap | F10's constants (G4) |
| material | `LadrunoSANISAND 1 264.32 0.312885 0.6944 1.3309 0.71 0.027 0.83 0.45 101 0.005 1.3 0.968 3.5 0.05 5.75 12.5 1100 2.0 <IntScheme> <TanType> 1 1e-7 <TolR> -flipAlphaIn init -Pmin 0.0101 -maxSubsteps 2000 -Presidual 0 -honorTolR 0` | intake §1 |
| staging | stage 0 (elastic) self-weight in 10 steps → stage 0 surface loads in 10 steps → `updateMaterialStage 1` → push | 86 emitter guide |
| threads | MKL = OMP = 1, `MKL_CBWR=COMPATIBLE`: bit-reproducible (A_r2 reproduced A bit-for-bit over 28 steps) and the same MKL path for both binaries | G5 |

Gravity controls (every leg): base reaction = γ'·A to 3e-16; 1-D K0 patch at the element average to 2.5e-12; η/M_c at the flip 0.769; min p' at the flip 5.39 kPa.

**Gaps filled (G1–G6).**
- **G1, mesh.** 22.5 × 18 m at a uniform B/8 would be 11 520 elements, not 2 430, so the spec contradicts itself. The act's own `ring_points_b8.csv` fixes it. Its element tags step by 12 per column along the surface over −5.7 ≤ x ≤ 3.1 m at 0.1875 m spacing, with tag 1950 at (0.8438, −0.0937). Only a 90-column mesh with a 15-row lower block numbered first fits that (90 × 27 = 2 430). With the reconstruction every tag in that CSV lands at its stated centroid. Not recoverable: where exactly the B/8 band ends in x (≥ 5.8 m; I took 6 m = 4B) and the grading law (geometric).
- **G2, unit weight.** §1 says "self-weight". The act's F10 intake (`f10_selfweight_wall.py` header) states γ' = 9.81 kN/m³ buoyant; used here.
- **G3, footing kinematics.** A rough, rigid, guided footing: u_x and θ of the reference node fixed. The alternative is a free rotation, which lets the bands pick a side. The act saw asymmetric bands; on this single-threaded, perfectly symmetric deck the fields stay symmetric to the last digit so far (the ring lists pair ±x points).
- **G4, controller.** The act's step floor and growth law are not in the spec; F10's (the fork's transcription of the same deck) are used.
- **G5, determinism.** Single-threaded MKL instead of 8 threads: the intake §1.6 run-to-run spread of 30 % comes from the threaded Pardiso, which also costs ~7 % of a step here. Binary A predates WP-132's `-deterministic`, so `MKL_NUM_THREADS=1` + `MKL_CBWR=COMPATIBLE` is the way to get determinism on both binaries.
- **G6, σ_zz.** The plane-strain SANISAND wrapper returns σ_xx, σ_yy, σ_xy only (`getStressZZ` = NaN). The replay files need all six components, so σ_zz is **recovered exactly** from the committed ψ and e (p_r = 0: e_c = e − ψ = e0 − λc (p/Pa)^ξ ⇒ p, σ_zz = 3p − σ_xx − σ_yy). It is checked at every checkpoint against the material's own `yieldDistance`: |f_read − f_recomputed| ≤ 1.5e-12 on all 9 720 points (quirks row added).

## 2. Control — DruckerPrager on the same deck

UW `DruckerPrager`, φ = 38°, ψ = 0 (ρ̄ = 0), plane-strain cone matched at ψ = 0 (√J2 = sin φ · p ⇒ ρ = √2 sin φ / 3), G = 30 MPa, K from ν*.
`updateMaterialStage` does not reach UW DP (existing quirks row); its K0 state is inside the cone, so it stays elastic anyway.

**Result: TARGET s/B 0.15 reached in 246 s**: 278 steps, 7 ladder failures absorbed by halving, rungs N/LS/K = 258/12/8. q peaks at 799 kPa at s/B 0.086, softens to ~710 kPa and recovers to 824 kPa at 0.145 (non-associated perfect plasticity). The act's §1.1 (every DP leg to 0.15 in 35–139 s, 0–3 subdivisions) is reproduced in kind, so the deck is sound. (`runs/ctrl_dp38/`)

## 3. Baseline — LadrunoSANISAND, ModifiedEuler (binary A, build 234a75751, WP-127 counters)

Leg `runs/A_me_baseline/`: IntScheme 1, TanType 0, TolR 1e-7 (inert under `-honorTolR 0`: TolE = 1e-4).

- **Reached s/B 0.01737, q 441.6 kPa, 79 steps, then stopped on the WALL-CLOCK cap (36 000 s), NOT on the step floor.** The machine paused ~2.75 h inside that budget: every run stalled together from 15:45 to 18:30. 6 failed attempts, all absorbed. 69 ModifiedEuler cap hits (2 000 substeps) on trial iterates, none on a converged step.
- The rerun (`A_me_baseline_r2`) matched it bit for bit for 28 steps. It was then stopped in the owner's re-plan, so **the act's wall (s/B 0.026–0.041) was not reached on this machine**. At s/B 0.0174 there was no sign of it: the controller was still at ds 2.5e-4–1e-3 m.
- **Ring at the last converged step:** the highest ρ_α (0.999) sits in the TOP ROW just outside the footing edge (elements 1842/1950, x = ±0.79 m, y = −0.04 m, p' 115 kPa). ρ_α 0.98–0.99 extends down the edge column to y = −0.34 m at p' 44–51 kPa. The lowest p' (3.97 kPa) is at the surface at x = ±2.5 m (x/B 1.7), with ρ_α 0.81–0.82. No point has ρ_α > 1.
- **Cost:** 9 720 points, ~3 000 substeps per point per step (summed over Newton iterations), 1–4·10⁷ substeps per step, 3–30 min per step at ds 1e-3 m. See §6 for the per-point split.

## 4. SAS-ME (binary B)

Builds, in order:
1. `cdf43685f`. Source-identical to WP-129 head 5c8dcd0e0: `git diff -- SRC` is empty. **Provisional**: it carries the known α-bounding dead-end and a non-error-controlled elastic path. Leg `runs/B_sasme_provisional_cdf43685f/`, stopped at step 45 (s/B 0.01104, q 311.9 kPa). **Dead-end census on it: zero.** No refusal of any kind, no committed ρ_α > 1.1 (max 0.966), no step cut. The defect never fired on this deck up to s/B 0.011.
2. **`beb6d8333`** (WP-129 review fixes, snapshot `binB_beb6d8333/`, md5 420e34d2…). **This is the result that counts.** Leg `runs/B_sasme_beb6d8333/`: IntScheme 129, TanType 0, TolR 1e-4, `-alphaEntryTol` 2 (default), κ 0.1, `-maxSubsteps` 2000.

**B, final local state.** B was stopped by the orchestrator (owner: moved to Esmeralda) after its last completed **step 70: s/B 0.017373, q 443.28 kPa**. That is the same s/B at which A ended, where A had 441.58 kPa.
- 4 ladder failures, absorbed.
- 134 refusals, all `maxSubsteps`, in ONE step (58). They are on 70 points in the TOP ROW at x = +1.8–2.0 m (x/B 1.2–1.35, the surface ring) and were absorbed by the ladder.
- 0 dead-ends; max ρ_α < 1.

**Cross-platform handover (E_B on Esmeralda, build 7936ed6e0, vs local B, beb6d8333).**
- Identical s/B sequence and q to 1e-5 kPa through step 40 (0.053 kPa by step 46).
- The step control forked at step 47: E_B converged on the Krylov rung at ds 1e-3 m; B failed one attempt and cut to 5e-4 m.
- Cost to step 65: E_B 0.98 h, B 5.59 h.

| s/B | E_B q (kPa) | B q (kPa) | E_B − B |
|---|---|---|---|
| 0.016040 | 416.787 | 416.801 | −0.014 |
| 0.016373 | 423.334 | 423.321 | +0.013 |
| 0.017040 | 436.436 | 436.708 | −0.27 (both real steps, after the fork) |

E_B reached step 65 at s/B 0.017373, q 443.145 kPa, with 0 refusals and 0 failed attempts after step 48. From here the WP-138 footing results come from the Esmeralda legs (E_A / E_B / E_C2 / E_B16).

## 5. Comparison (so far)

### 5.1 q–s, B vs A (`compare_curves.py A_me_baseline B_sasme_beb6d8333`)

| s/B range | max \|q_B − q_A\| / q_A |
|---|---|
| 0 – 0.001 | **5.5 %** (at the first step: A accepts s/B 1.3e-5 on the LineSearch rung at 13.76 kPa, B at 13.00; both are ~12–14 kPa) |
| 0.001 – 0.005 | 0.64 % |
| 0.005 – 0.0174 | 0.51 % (worst absolute 2.05 kPa at s/B 0.0170; at the common end 443.28 vs 441.58 kPa, +0.38 %) |

An earlier message quoted "within 0.31 %". That figure is the worst **absolute** difference (1.1 kPa) over s/B ≤ 0.0134, taken as relative at that point. It is not the worst relative difference, and it does not hold below s/B 0.005. The table above is the corrected statement.

Plots: `Ladruno_files/testbed/footing_sas_me_ab/qs_local_zoom.png` (A, B, DP to s/B 0.018) and `qs_local.png` (all local legs). The wall-clock panel has jumps, and they are not integrator cost:
- A's jump near s/B 0.009 is the 2.75 h machine pause.
- B's jump near s/B 0.012 is one 49-min step, 47, which ran while four runs shared the machine.

### 5.2 Cost to the common end s/B 0.01737

| leg | push wall | substeps | Newton iterations | steps | failed attempts | caps / refusals |
|---|---|---|---|---|---|---|
| A (ModifiedEuler) | 10.05 h (≈ 7.3 h without the 2.75 h pause) | 1.22e9 | 1 514 | 79 | 6 | 148 cap hits (trial iterates) |
| B (SAS-ME beb6d8333) | 7.03 h | 1.24e9 | 1 287 | 70 | 4 | 134 refusals (one step) |

B's wall clock also shares the machine with up to three other runs; the substep and iteration counts are the load-independent measure.

On this deck SAS-ME costs about the same substeps per unit settlement as ModifiedEuler. The expected 1.5–3× substep overhead per increment did not materialise at the global level.

### 5.3 Material-point replays with the REAL increments (reconciled)

Each replay row is one Gauss point's COMMITTED state from the run plus the strain increment that point actually received in the next converged step (`dStrain`, compression positive, engineering γ). Every dumped pair is replayed by the binary that produced the run to ≤ 3e-12 in σ, so the dumps are faithful. Rows are replayed through:
- ModifiedEuler (binary A, build 234a75751);
- SAS-ME (binary B, beb6d8333 unless stated);
- the WP-134 independent oracle (`sanisand_reference`, preset `uw_model`, Radau rtol 1e-10).

**Error norm used everywhere:** e = ‖σ_C++ − σ_oracle‖ / p'_in (contravariant norm, all 6 components). Produced by `run_replays.py` + `replay_summary.py`; inputs in `runs/<leg>/replay/`, outputs in `runs/<leg>/replay_out/`.

| committed state (A run) | Δ of the next step | n | ME median / p90 / max | SAS-ME median / p90 / max | refusals |
|---|---|---|---|---|---|
| step 25, s/B 0.00141 | ds 3.2e-4 m | 44 | 1.16e-3 / 3.53e-3 / 4.30e-3 | 5.74e-5 / 7.58e-5 / 1.10e-4 | 0 / 0 |
| step 50, s/B 0.01246 | ds 1.25e-4 m | 46 | 6.37e-5 / 9.87e-5 / 3.22e-4 | 6.28e-5 / 7.79e-5 / 1.12e-4 | 0 / 0 |
| step 75, s/B 0.01654 | ds 5.0e-4 m | 43 | 9.42e-5 / 1.30e-4 / 1.43e-4 | 1.85e-4 / 2.01e-4 / 2.17e-4 | 0 / 0 |
| step 78 → 79 (last converged pair), s/B 0.01721 | ds 2.5e-4 m | 44 | 6.85e-5 / 8.55e-5 / 1.98e-3 | 9.99e-5 / 1.15e-4 / 1.25e-4 | 0 / 0 |
| step 25 of the cdf43685f run (SAS-ME by cdf43685f) | ds 3.2e-4 m | 44 | 1.14e-3 / 3.46e-3 / 4.21e-3 | 7.22e-5 / 8.85e-5 / 1.08e-4 | 0 / 0 |

**Reconciliation.**
- The earlier "ME median 1.2e-3·p', SAS-ME 6e-5·p'" was the **step-25** dataset: early loading, s/B 0.0014, where the ring points are first becoming plastic. It is committed, as `replay_step00025.*`.
- The drafter's "~7e-5 / ~1e-4" is the **last-converged-pair** dataset at s/B 0.017.
- Both are right for their own dataset. Same binaries, oracle preset and norm; the SAS-ME step-25 row was re-run on beb6d8333 and agrees with the cdf43685f figure.

**The accuracy claim that survives.**
- **SAS-ME** sits at 0.6–2e-4·p' at every checkpoint: its TolR 1e-4 error control, on the substeps.
- **ModifiedEuler**'s error depends on the state. It is 20× worse than SAS-ME at the onset of ring plasticity (step 25). It is equal or slightly better once the ring is fully mobilised (steps 50–78), except for isolated low-p outliers (2e-3·p' at the last pair, a point taken in one substep).
- **A blanket "20× more accurate" is WRONG.** State it as "20× at the onset of plasticity (s/B ≈ 0.001), parity from s/B ≈ 0.01".
- **No real increment was refused** by either integrator, and the oracle integrated all of them (`ok`). Up to s/B 0.017 the points' actual increments are integrable.

## 6. Where the cost goes — per-point substep split (`census_dist.py`)

Ring definitions: **ringP** = p' < 10 kPa (830–850 points); **ringG** = top element row with 0.5 ≤ |x|/B ≤ 2 (96 points).

| leg | interval (s/B) | substeps / step | ringP share | ringG share | top-1 | top-10 | top-100 |
|---|---|---|---|---|---|---|---|
| A | 0.0014–0.0125 | 3.0e7 | 5.4 % | 3.1 % | 0.10 % | 0.98 % | 7.4 % |
| A | 0.0125–0.0165 | 1.5e7 | 5.4 % | 2.5 % | 0.13 % | 1.3 % | 8.9 % |
| A | 0.0165–0.0174 | 2.3e7 | 6.2 % | 2.9 % | 0.13 % | 1.3 % | 8.7 % |
| B | 0.0014–0.0129 | 3.1e7 | 7.7 % | 3.1 % | 0.09 % | 0.87 % | 6.8 % |

- Per-step top-1 share over all steps: A median 0.14 %, max 0.34 %; B median 0.12 %, max 0.32 %.
- Histograms: ~99 % of the points take 1e3–1e5 substeps per interval.
- **The cost is spread over the whole mesh, not the ring.** A per-point SAS-ME → CPPM fallback for the worst points would recover < 10 %. The multiplier is the Newton iteration count × every point's update, so the lever is fewer global iterations (the TanType 1 legs C2/D, run elsewhere) or a threaded update (F19).

## 7. Oracle-replayable checkpoints (format)

`runs/<leg>/replay/replay_stepNNNNN.csv`:
- The ~50 worst points at a checkpoint: top 20 by ρ_α, the 15 lowest p', and the top 15 by substeps in the next step.
- Columns: the TIMs ring CSV's columns (element, gp, x_m, y_m = element centroid, p_kPa, eta, eta_over_Mb_compression, e, psi, sigma_0..5, alpha_0..5, alpha_in_0..5, z_0..5), then `dStrain_0..5`, step, s_over_B, gp_x_m, gp_y_m, rho_alpha, f_read, f_recomputed, substeps_next, capHit_next, dt_next, prevIncrNorm, sigma_next_0..5, alpha_next_0..5, e_next, select.
- **Sign:** σ, dStrain and σ_next are COMPRESSION positive (feeds `ladrunoSANISANDReplay -convention compressionPositive` and `sanisand_reference` unchanged). α, α_in and z are the raw internal ratios.
- `replay_wall_last_pair.csv` = the state at step n−1 plus the increment of the last converged step n.
- `replay_wall_probe_iter1.csv` (written only at a FLOOR wall) = the state at n plus the first Newton iterate of the failing increment, committed through a `FixedNumIter 1` post-mortem (quirks row).

## 8. Verified / not verified

- **Verified:**
  - the gravity patch and resultants;
  - σ_zz recovery (1.5e-12);
  - self-replay of every dump by its own binary (≤ 3e-12);
  - determinism (A_r2 = A bit for bit);
  - the DP control reaches 0.15.
- **Not verified:**
  - the act's wall location (not reached: see §3);
  - the TanType 1 legs (C2/D, another agent);
  - SAS-ME to its end: pending.

## 9. Pending

B's end state; C2/D (TanType 1, Krylov at 1× tolerance, TolR 1e-4 / 1e-3) from the orchestrator; the final q–s plot (`plot_qs.py`); the verdict.
