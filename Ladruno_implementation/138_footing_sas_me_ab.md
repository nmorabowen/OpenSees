---
title: "WP-138 — strip-footing A/B: ModifiedEuler (IntScheme 1) vs SAS-ME (IntScheme 129)"
project: Ladruno
type: measurement report
status: MEASURED — the Esmeralda arms are final (§8–§10, verdict §0); the sensitivity ladders are an INTERIM snapshot (§11); the default in §12 is a recommendation for the owner/TIMs to decide
related:
  - "[[_tims_2d_model_requests_2026-09-25]]"
  - "[[LadrunoSANISAND_implex_guide]]"
  - "[[150_sanisand_regularization_memo]]"
  - "[[LEDGER_quirks]]"
updated: 2026-09-28
---

# WP-138 — strip-footing A/B, ModifiedEuler vs SAS-ME (TIMs F18 end to end)

Pure Python decks and runs; no C++ change, no build. Deck, driver and every run
output: `Ladruno_files/testbed/footing_sas_me_ab/`. Nothing in the TIMs Workbench
was run or edited; the deck is built from the intake's §1 spec
(`_tims_2d_model_requests_2026-09-25.md`, on `origin/wp/127-sanisand-replay-counters`).

## 0. Verdict

**SAS-ME moves the wall from s/B 0.0292 to 0.0508, but it does not remove it. The
wall is a constitutive singularity, not an integration defect. There is no peak and
no plateau on this deck.**

- **E_A, ModifiedEuler** (IntScheme 1): stops on the step floor at **s/B 0.0292, q 701.8 kPa**. That is inside TIMs' own ModifiedEuler wall band (s/B 0.026–0.041), so the fork reproduces TIMs. On the way in, ModifiedEuler force-accepts at dt_min and commits ρ_α up to 13.09. The result is a spurious +6 % stiffening against SAS-ME at the same s/B (§8.4).
- **E_B, SAS-ME** (IntScheme 129, TanType 0, TolR 1e-4): stops on the step floor at **s/B 0.0508, q 966.7 kPa**. q is still rising there: the slope over the last 0.005 s/B is 0.24 × the initial slope, and q_max = q_end. The first `loadingNonPosH` refusal comes at s/B 0.0363. From then on NonPosH refusals accumulate, and they end the run.
- **Every SANISAND arm stops in MODE = FLOOR.**
  - The SAS-ME arms (E_B, E_D, E_C2, E_B16) reach the floor on `loadingNonPosH` refusals.
  - ModifiedEuler has no refusal path, so E_A reaches it through cap hits and forced-at-dt_min acceptances at the footing-edge point.
  - The DruckerPrager control on the same deck runs to s/B 0.15 (§2).
- **The mechanism** (WP-150 memo, draft #892; §10) is DM04's hardening-modulus singularity at α_in re-seats. At a = (α − α_in):n → 0, h reaches its 1e10 cap. With b:n ≤ 0 this sends K_p → −∞.
  - The refusing points are PRE-peak (ρ_α < 1) and they chatter.
  - The WP-134 oracle hits the same 0/0.
  - No integrator knob lifts it: TolR 1e-3 walls EARLIER (0.0410), and the consistent tangent walls at 0.0114.
  - Neither does a BVP regularizer.
- **Default** (§12; recommendation, owner/TIMs decide): IntScheme 129, TanType 0, TolR 1e-4 with the step policy these runs used. **Follow-ups** (§13): R1 is a CONSTITUTIVE change (bounded h), and that decision belongs to the owner/TIMs.

| arm | integrator | s/B at FLOOR | q (kPa) | first NonPosH s/B | refusals (converged-step census) | push wall (h) |
|---|---|---|---|---|---|---|
| E_A | ModifiedEuler | **0.0292** | 701.8 | — (no refusal path) | 542 cap hits, 3 943 forced at dt_min | 2.90 |
| E_B | SAS-ME, TolR 1e-4 | **0.0508** | 966.7 | **0.0363** | NonPosH 232, maxSubsteps 209, errorAtDTmin 1 | 4.32 |
| E_D | SAS-ME, TolR 1e-3 | 0.0410 | 808.3 | 0.0334 | maxSubsteps 13 316, errorAtDTmin 190, NonPosH 173 | 4.67 |
| E_C2 | SAS-ME, TanType 1 | 0.0114 | 317.0 | 0.0064 | errorAtDTmin 175 011, maxSubsteps 40 683, NonPosH 409 | 12.35 |
| E_B16 | SAS-ME, B/16 mesh | 0.0135 | 352.8 | 0.0135 | maxSubsteps 115, NonPosH 11 | 3.75 |
| ctrl DP 38° | DruckerPrager | 0.15 (TARGET) | 752.0 (max 824.2) | — | — | 0.07 |

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
- **Superseded for the wall by E_A on Esmeralda (§8):** the same integrator on build 7936ed6e0 walls at s/B 0.0292, inside the act's band.

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

## 5. Comparison — local legs, s/B ≤ 0.0174

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

**Reconciliation — which increments each figure covers.**
- **"ModifiedEuler 1.2e-3·p′ vs SAS-ME 6e-5·p′"** are the medians of the **step-25** set, `runs/A_me_baseline/replay_out/replay_step00025.*`.
  - It uses the A run's committed state at step 25 (s/B 0.00141).
  - It covers 44 selected points (top ρ_α, lowest p′, top substeps).
  - Each point is replayed through the real increment it received in step 26 (ds 3.2e-4 m).
  - ME 1.16e-3; SAS-ME 5.74e-5 on beb6d8333 (7.22e-5 on cdf43685f, from that run's own step 25).
  - This is early loading, where the ring points first become plastic.
- **"≈7e-5 / ≈1e-4"** are the medians of the **last-converged-pair** set, `replay_wall_last_pair.*`.
  - It uses the state at step 78 plus the increment of step 79 (s/B 0.01721, ds 2.5e-4 m).
  - It covers 44 points.
  - ME 6.85e-5, SAS-ME 9.99e-5. So at that pair the ModifiedEuler **median** is the lower one; its max, 1.98e-3, is not.
- Both figures are right for their own increments. They share the binaries, the oracle preset and the norm. They are medians over ~44 selected worst points, **not** over the 9 720 points of the mesh.
- In between: at step 50 (s/B 0.0125) the two are at parity, 6.37e-5 vs 6.28e-5. At step 75 (s/B 0.0165) ME is lower, 9.42e-5 vs 1.85e-4.
- **No replay exists at the Esmeralda walls.** The E_* runs wrote `replay_wall_last_pair.csv` on Esmeralda, but none was replayed. Every replay figure therefore covers s/B ≤ 0.0174.

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
- **The cost is spread over the whole mesh, not the ring.** A per-point SAS-ME → CPPM fallback for the worst points would recover < 10 %. The multiplier is the Newton iteration count × every point's update, so the lever is fewer global iterations or a threaded update (F19).
- **Fewer global iterations did not come from the consistent tangent.** E_C2 (TanType 1, §8.6) cut the iterations per step, but its steps collapsed to 1e-5 to 1.25e-6 m. It cost 13× E_B per unit s/B.

## 7. Oracle-replayable checkpoints (format)

`runs/<leg>/replay/replay_stepNNNNN.csv`:
- The ~50 worst points at a checkpoint: top 20 by ρ_α, the 15 lowest p', and the top 15 by substeps in the next step.
- Columns: the TIMs ring CSV's columns (element, gp, x_m, y_m = element centroid, p_kPa, eta, eta_over_Mb_compression, e, psi, sigma_0..5, alpha_0..5, alpha_in_0..5, z_0..5), then `dStrain_0..5`, step, s_over_B, gp_x_m, gp_y_m, rho_alpha, f_read, f_recomputed, substeps_next, capHit_next, dt_next, prevIncrNorm, sigma_next_0..5, alpha_next_0..5, e_next, select.
- **Sign:** σ, dStrain and σ_next are COMPRESSION positive (feeds `ladrunoSANISANDReplay -convention compressionPositive` and `sanisand_reference` unchanged). α, α_in and z are the raw internal ratios.
- `replay_wall_last_pair.csv` = the state at step n−1 plus the increment of the last converged step n.
- `replay_wall_probe_iter1.csv` (written only at a FLOOR wall) = the state at n plus the first Newton iterate of the failing increment, committed through a `FixedNumIter 1` post-mortem (quirks row).

## 8. The Esmeralda arms (final)

**Where the runs and their records are.**
- Build: ladruno **7936ed6e0** (`esmeralda/build_wp138.sh`, Linux oneMKL 2024.2).
- Deck: this deck, `esmeralda/footing_ab_esmeralda.py`, with `--mesh b16` for E_B16.
- Each arm ran alone on its node, MKL/OMP = 1 thread. Wall clock is therefore comparable between arms (not with the local legs).
- Launch record: `esmeralda/JOBS.txt`.
- Run records: `runs/E_A`, `E_B`, `E_D`, `E_C2`, `E_B16`, each with `steps.csv`, `summary.json`, `logs/log.log` and `rung_fail.csv`. E_C2 was re-pulled after its end (exit 0 after 44 452 s).
- Tables: `python wall_table.py`.
- Figures: `python plot_qs_final.py`.
- The census/field analyses: `esmeralda/analysis/` (`a1_walls.py`, `a2_curves_fields.py`, `a3_refusers.py`, `tables/`, PNGs). They were run at 12:43, before E_C2 ended, so their E_C2 entries and the "(running)" labels are stale. For E_C2 use `runs/E_C2` and `qs_final*.png`.
  - They read the checkpoint `.npz` fields, which stay on Esmeralda (`~/ladruno_wp138/deck/runs/<leg>/ckpt/`) and are not committed.
  - The same goes for `census_last_converged.csv`, 1–7 MB per arm.

| arm | settings (differ from E_B) | mode | s/B end | q end (kPa) | steps | failed attempts | first NonPosH s/B (step) | refusals by code, converged-step census | push wall (h) | h per 0.01 s/B |
|---|---|---|---|---|---|---|---|---|---|---|
| E_A | IntScheme 1, TolR 1e-7 (inert) | FLOOR | 0.02923 | 701.8 | 189 | 31 | — | ME: 542 cap hits / 3 943 forced at dt_min | 2.90 | 0.99 |
| E_B | IntScheme 129, TanType 0, TolR 1e-4 | FLOOR | 0.05084 | 966.7 | 377 | 57 | 0.03629 (180) | loadingNonPosH 232, maxSubsteps 209, errorAtDTmin 1 | 4.32 | 0.85 |
| E_D | TolR 1e-3 | FLOOR | 0.04100 | 808.3 | 389 | 57 | 0.03342 (185) | maxSubsteps 13 316, errorAtDTmin 190, loadingNonPosH 173, tensionAtDTmin 1 | 4.67 | 1.14 |
| E_C2 | TanType 1, -maxSubsteps 20 000, Krylov tol ×1 | FLOOR | 0.01144 | 317.0 | 2 283 | 345 | 0.00641 (342) | errorAtDTmin 175 011, maxSubsteps 40 683, loadingNonPosH 409 | 12.35 | 10.77 |
| E_B16 | B/16 mesh (38 880 GP) | FLOOR | 0.01352 | 352.8 | 77 | 17 | 0.01350 (76) | maxSubsteps 115, loadingNonPosH 11 | 3.75 | 2.77 |
| ctrl_dp38 | DruckerPrager 38°, local Windows run | TARGET | 0.15 | 752.0 (max 824.2) | 278 | 7 | — | — | 0.07 | 0.005 |

**How the table was counted.**
- The census sums the `refusals step N: … by code {…}` lines. They are written for every converged step and include that step's failed attempts. The final ladder at the floor never converges, so it is not in the sum.
- The floor-ladder census is in `esmeralda/analysis/tables/walls_summary.json`:
  - E_B: +26 NonPosH at 2 points.
  - E_B16: +24 NonPosH at 1 point.
  - E_D: +24.
  - E_A: 24 cap hits.
- "h per 0.01 s/B" = wall clock at the last converged step ÷ (s/B_end / 0.01).

### 8.1 q–s to the end of every arm

`Ladruno_files/testbed/footing_sas_me_ab/qs_final.png` (full) and `qs_final_zoom.png` (s/B ≤ 0.06). In the plot, ○ marks where each arm stopped and ▽ marks the first NonPosH step.
- **E_B has no peak and no plateau.**
  - q_max = q_end = 966.7 kPa.
  - The end slope is 0.24 × the initial slope over the last 0.005 s/B, and 0.31 × over the last 0.001.
  - One 4.35 kPa dip (0.5 %) at s/B 0.0490–0.0500 recovers by 0.0500. After it q rises another 12 kPa to the wall.
- **E_A turns UP before its wall.** Its end slope is 1.37 × the initial slope over the last 0.001 s/B. At its wall it is **+6.4 %** above E_B at the same s/B (§8.4).
- **The ModifiedEuler and SAS-ME curves agree until the ME wall approaches.**
  - The Esmeralda pair reproduces the local A-vs-B figures (§5.1): max |Δq|/q is 0.65 % over s/B 0.001–0.005 and 0.54 % over 0.005–0.016.
  - Over 0.016–0.029 the gap grows to 6.1 %. The growth is E_A turning up while it commits states outside the bounding surface (§8.4).
- **The DP control is not a SANISAND bound.** It is a matched-cone plasticity check that the deck can carry a mechanism to 0.15 (§2).

### 8.2 Wall clock and substeps per unit s/B

Hours / 1e9 substeps per 0.01 s/B, by interval (`wall_table.py`; a partial interval is prorated):

| arm | 0–0.01 | 0.01–0.02 | 0.02–0.03 | 0.03–0.04 | 0.04–0.05 |
|---|---|---|---|---|---|
| E_A (ME) | 0.58 / 0.48 | 1.14 / 0.94 | 1.24 / 1.01 | — | — |
| E_B (SAS-ME) | 0.40 / 0.49 | 0.82 / 1.02 | 0.93 / 1.09 | 1.02 / 1.14 | 1.05 / 1.14 |
| E_D (TolR 1e-3) | 0.47 / 0.44 | 0.99 / 0.97 | 0.98 / 0.96 | 1.43 / 1.39 | 7.82 / 7.88 (0.040–0.041 only) |
| E_C2 (TanType 1) | 5.12 / 5.05 | 50.0 / 56.6 (0.010–0.0114 only) | — | — | — |
| E_B16 (B/16) | 2.46 / 2.71 | 3.65 / 4.16 (0.010–0.0135 only) | — | — | — |

- **Substeps per unit settlement: SAS-ME ≈ ModifiedEuler** (0.48 vs 0.49 and 0.94 vs 1.02 per 0.01 s/B). This confirms the local finding (§5.2).
- **Wall clock: SAS-ME is 1.3–1.45× cheaper per unit s/B** (0.40 vs 0.58, 0.82 vs 1.14, 0.93 vs 1.24 h).
- **E_B's cost is flat from s/B 0.01 to its wall**, about 1 h per 0.01 s/B, so the 0.0508 wall is not a budget stop.
- **Krylov carries the settlement.** The KrylovNewton rung at tol × 10 carries 80–90 % of the settlement in every arm: `frac_settlement_on_K` in `walls_summary.json` is E_A 0.89, E_B 0.91, E_D 0.93, E_B16 0.82.
  - Accepting at 10× the tolerance costs about ±1.5 kPa on q.
  - The deck's own curve gate is 0.028 kPa (0.0415 kN / 1.5 m).
  - Read the curves with that band.

### 8.3 E_B — SAS-ME to its wall

- The first `loadingNonPosH` comes at step 180, s/B 0.03629 (element 1829, gp 3, at (−0.98, −0.23)). The refusals keep coming from there.
- The floor ladder refuses at **2 points under the footing core**:
  - element 1880 gp 1 at (−0.148, −2.02) m, p′ 366 kPa, **ρ_α 0.958**, 27 refusals;
  - element 1879 gp 1 at (−0.148, −2.21) m, 3 refusals.
- Both are **pre-peak**, and both sit in the band (`floor_refusers_in_band.csv`: 98th percentile of total γ).
- At the last converged step: 81 points with ρ_α > 1, max 1.007; p′_min 2.17 kPa; none below 1 kPa.
- The failed rungs are 868 in all:
  - Newton divergence 391;
  - LineSearch divergence 281;
  - **LineSearch refusal 119**;
  - Krylov divergence 36;
  - Newton refusal 20;
  - Krylov refusal 21.

### 8.4 E_A — ModifiedEuler at TIMs' wall

- **The floor.** The arm floors at s/B 0.02923 on the footing-edge point: element 1950 gp 2 at (0.898, −0.148) m, p′ 58.9 kPa, ρ_α 1.0008. At the floor that point shows 31 cap hits and 513 forced-at-dt_min acceptances.
- **No refusal path, so bad states get committed.** ModifiedEuler has no refusal path, so it commits what it cannot integrate. Over s/B 0.026–0.0293 it forced 3 786 acceptances at dt_min, and **the committed ρ_α reaches 13.09**. SAS-ME, over the same window: max 1.004 and no forced acceptances.
- **The curve turns up.** The curve bends UP to +6.4 % over E_B at the wall. That stiffening is spurious; it comes from the accepted states outside the bounding surface.
- **Reading TIMs' curves.** TIMs' ModifiedEuler curves near their wall should be read with this in mind.

### 8.5 E_D — TolR 1e-3 is not the lever

- **E_D walls EARLIER than E_B:** 0.0410 against 0.0508. Keep TolR 1e-4.
- **From the E_D records:**
  - Its cost per unit s/B equals E_B's up to s/B 0.038: 0.98 vs 0.93 h per 0.01 in 0.02–0.03.
  - Near its wall the cost rises to 4.6× E_B's (0.038–0.041) and 8× (0.040–0.041).
  - It had 13 316 maxSubsteps refusals against E_B's 209.
  - Its q runs **−3.5 % (median) below E_B** over s/B 0.02–0.041, with a range of −4.3 % to +1.6 %.
- **The campaign's own TolR study,** from the orchestrator (not in this branch's records), reports TolR 1e-3 as **27× slower and +0.54 % biased**, and refutes the "2× faster" lever.
  - The two sets of figures measure different things. Their conclusion is the same: **no saving, and an earlier wall.**
  - The source of the 27× / +0.54 % figures should be named wherever they are quoted.

### 8.6 E_C2 — the consistent tangent (TanType 1) is not viable here

E_C2 is E_B with TanType 1, `-maxSubsteps 20000` and the KrylovNewton rung at 1× tolerance.

**How the steps collapsed:**
- Newton iterations per step fell: median 11 → 4 → 2.
- The accepted step size collapsed with them: median ds 2e-5 m to s/B 0.005, 1e-5 m to 0.01, 1.25e-6 m after that. E_B runs at 1e-3 m there.
- 2 283 steps; 175 011 errorAtDTmin and 40 683 maxSubsteps refusals.
- The first NonPosH came at s/B 0.0064, the earliest of any arm.
- It floored at s/B 0.0114 after 12.35 h: 10.8 h per 0.01 s/B, **13× E_B**, and 61× in 0.01–0.0114.

**Its curve still tracks E_B** within 1.45 % to its end. The consistent tangent does not buy settlement; it buys a step-size collapse. Keep TanType 0. The two changes were made together, so the cost of Krylov at 1× on its own is **not measured**.

## 9. B/16 — non-associated localization, mesh-dependent

E_B16 is E_B on a B/16 mesh: 180 × 54 = 9 720 elements (38 880 Gauss points). It has the same domain, fine-band extents, loads and BCs; every B/8 cell becomes ~2 × 2 cells, and there are 17 footprint nodes. It floors at **s/B 0.0135**, q 352.8 kPa. Its first NonPosH comes at 0.0135 (step 76), and the floor ladder refuses at one point (element 7820 gp 4, at (0.957, −0.395) m, p′ 50 kPa, ρ_α 0.93, pre-peak).

**The curve.** B/16 runs softer than B/8 from s/B 0.010 (`curves_fields_summary.json`):

| s/B | 0.002 | 0.005 | 0.008 | 0.010 | 0.012 | 0.013 | 0.0135 |
|---|---|---|---|---|---|---|---|
| q_B16 / q_B8 − 1 | −0.58 % | −0.14 % | −0.91 % | −3.89 % | −4.95 % | −4.84 % | −3.57 % |

**The band is one element wide on both meshes.**
- The FWHM of the incremental shear strain across the band at four depths is 1.01–1.11 h on B/16 (0.094–0.104 m) and 1.00–1.05 h on B/8 (0.19–0.20 m).
- The band halves with the element, it follows the mesh lines, and the total-γ peak sits at the footing edge (x ≈ ±0.8 m). Plot: `esmeralda/analysis/band_profiles_B8_vs_B16.png`.

**What it is.** It is **non-associated localization** (Rudnicki–Rice), **not ψ-softening**. Per #892 (the WP-150 memo):
- The plane-strain acoustic-tensor scan of the continuum tangent on the E_B/E_B16 checkpoints finds det ≤ 0 at 17 % of the Gauss points already at s/B 0.011 (B/8) and 0.0096 (B/16).
- At that point the material is still hardening (H/2G ≈ 1.05).
- The associated control (R → Q) is elliptic everywhere.
- Only 17 of 9 720 points are post-peak at s/B 0.0508.

That explains both the 4–5 % softening of B/16 from s/B 0.010 and the mesh-line bands. **The wall is a separate matter** (§10): the floor refusers are pre-peak on both meshes. The wall's s/B on B/16 is earlier than on B/8 (0.0135 vs 0.0508). The mesh dependence of the wall is not settled; R2 (§13) re-measures it.

## 10. Why the wall — the diagnosis (per #892, the WP-150 memo)

The wall is DM04's hardening-modulus singularity at an α_in re-seat.

**The singularity.**
- With a = (α − α_in):n, h = b0 / a. As a → 0, h hits its 1e10 cap.
- With b:n ≤ 0, K_p → −∞, and the update has no solution. `loadingNonPosH` is the SAS-ME refusal that names it.
- The singular set is **{a = 0, b:n ≤ 0}**.

**How the load path reaches it.** Dilation raises ψ. That contracts the bounding surface onto α, so b:n → 0⁻.
- With A0 = 0.001 (S4, §11) the bounding surface stays an attractor, b:n stays > 0, and the refusal does not occur.
- With h0 × 3 α reaches the bounding surface sooner, so the arm walls earlier (onset s/B 0.0091).
- Presidual and e_init leave the set intact (§11).

**The refusers.**
- They are PRE-peak (ρ_α < 1 at every floor refuser in §8.3 and §9) and they chatter: per #892, E_B makes 10.1 M re-seats.
- The WP-134 independent oracle hits the same 0/0.
- No integrator setting lifts it: TolR 1e-3 walls earlier (§8.5), TanType 1 walls at 0.0114 (§8.6), and ModifiedEuler walls earlier still and commits ρ_α 13 states (§8.4).
- No BVP regularizer lifts it either. A nonlocal ψ̄ or a crack band does not act on this onset, which is set at the material point by a = 0 with b:n ≤ 0.
- So it is **constitutive, not an integration defect**. Any change must be made in the model (R1, §13), and that is an owner/TIMs decision.

## 11. Sensitivity ladders — interim snapshot 2026-09-28 16:20

> **INTERIM.** Taken at 16:20. At that time S1–S4, L_pres_1/2/5/10/20, L_e_0p65 and L_e_0p80 were still RUNNING.
>
> **To refresh:** `cd Ladruno_files/testbed/footing_sas_me_ab/ladders_interim && python collect.py`. It pulls from Esmeralda and rewrites `ladder_table.txt` and `q_s_{presidual,einit,ablation,A0_h0}.png`. Then replace the table and the verdict paragraph below.
>
> **Driver:** `esmeralda/footing_ab_esmeralda.py` + `patch_driver.py` (`--presidual`, `--einit`) + `patch_driver2.py` (`--zmax --nb --nd --A0 --h0`). Launch record: `esmeralda/JOBS.txt`.
>
> **Setup:** every leg is E_B with one knob changed. The ablation S1→S4 is CUMULATIVE. The S and A0/h0 legs run 4 to a node (4 CPUs, 6 GB), so no wall clock is quoted. **The ladder's p′/η/ρ diagnostics are wrong on the Presidual legs and are not quoted.**

| ladder | leg | knob | status | s/B reached | q (kPa) | first NonPosH s/B | refusals: NonPosH / maxSubsteps (other) |
|---|---|---|---|---|---|---|---|
| — | E_B | reference (Presidual 0, e 0.6944, A0 0.05, h0 1.3) | FLOOR | 0.0508 | 966.7 | 0.0363 | 232 / 209 (errorAtDTmin 1) |
| Presidual | L_pres_0p5 | 0.5 kPa | FLOOR | 0.0303 | 682.9 | **0.0182** | 38 / 484 |
| Presidual | L_pres_1 | 1 kPa | running | 0.0400 | 836.0 | 0.0333 | 201 / 2 576 |
| Presidual | L_pres_2 | 2 kPa | running | 0.0427 | 886.6 | 0.0346 | 128 / 2 023 |
| Presidual | L_pres_5 | 5 kPa | running | 0.0519 | 1 062.8 | 0.0416 | 22 / 585 |
| Presidual | L_pres_10 | 10 kPa | running | 0.0590 | 1 243.5 | 0.0373 | 8 / 828 |
| Presidual | L_pres_20 | 20 kPa | running | 0.0711 | 1 589.7 | 0.0535 | 21 / 1 |
| e_init | L_e_0p65 | 0.65 | running | 0.0435 | 1 391.9 | 0.0395 | 15 / 606 |
| e_init | L_e_0p75 | 0.75 | FLOOR | 0.0431 | 500.1 | 0.0349 | 99 / 1 412 |
| e_init | L_e_0p80 | 0.80 | running | 0.0567 | 337.8 | 0.0525 | 17 / 1 760 |
| e_init | L_e_0p85 | 0.85 | FLOOR | 0.0442 | 180.7 | 0.0414 | 143 / 1 040 (errorAtDTmin 5) |
| ablation | S1_nofabric | z_max = 0 | running | 0.0368 | 779.4 | 0.0355 | 8 / 948 |
| ablation | S2_nopeak | + n_b = 0 | running | 0.0348 | 316.1 | 0.0269 | 22 / 517 |
| ablation | S3_critstate | + n_d = 0 | running | 0.0362 | 303.0 | 0.0347 | **1** / 1 037 |
| ablation | S4_nodilat | + A0 = 0.001 | running | 0.0374 | 311.4 | **none** | **0** / 1 276 (errorAtDTmin 2) |
| A0 / h0 | A0_0p02 | A0 = 0.02 | FLOOR | 0.0315 | 677.7 | 0.0237 | 86 / 429 |
| A0 / h0 | A0_0p10 | A0 = 0.10 | FLOOR | 0.0361 | 834.8 | 0.0310 | 55 / 67 (errorAtDTmin 1) |
| A0 / h0 | h0_x3 | h0 = 3.9 | FLOOR | 0.0188 | 769.6 | **0.0091** | 222 / 3 257 |

The refusal counts come from two sources:
- FLOOR legs (other than E_B): the driver's `summary.json` census.
- RUNNING legs and E_B: the converged-step log census.

**Interim verdict (16:20).**
- **Only killing the dilatancy (S4) clears `loadingNonPosH`.**
  - It persists through S1 (fabric off), S2 (+ no peak) and S3 (+ critical-state dilatancy: 1 event).
  - It is ABSENT in S4 (+ A0 = 0.001) at s/B 0.0374, which is past E_B's onset (0.0363).
  - S4 is not cheap. Per the campaign it costs ~27 M substeps per step (the committed interim records give a median of 2.4e7 over its last 20 steps, against 7.3e6 for E_B), with maxSubsteps refusals. That cost is the re-seat chatter (R1b, §13).
- **Presidual 0.5–20 kPa never clears NonPosH.**
  - Its onset is non-monotonic: 0.5 kPa brings it to 0.0182, EARLIER than Presidual 0.
  - Larger Presidual walls later and stiffens q (q at s/B 0.03: 675 → 835 kPa from 0.5 to 20 kPa). That is an apparent cohesion, not a cure.
- **Every e_init from 0.65 to 0.85 shows NonPosH.**
- **A0 is non-monotonic:** onset at 0.0237 (A0 0.02), 0.0363 (0.05) and 0.0310 (0.10).
- **h0 × 3 brings the onset down to 0.0091.**
- This matches the §10 mechanism: the set {a = 0, b:n ≤ 0} is left intact by every knob except the one that keeps b:n > 0.

## 12. Default integrator for the TIMs campaign — recommendation; owner/TIMs decide

**Recommendation: IntScheme 129 (SAS-ME) with TanType 0 and TolR 1e-4, on the step policy these runs used.** As the E_B material line:

```
LadrunoSANISAND $tag 264.32 0.312885 0.6944 1.3309 0.71 0.027 0.83 0.45 101 0.005 1.3 0.968 3.5 0.05 5.75 12.5 1100 2.0 \
    129 0 1 1e-7 1e-4 -flipAlphaIn init -Pmin 0.0101 -maxSubsteps 2000 -Presidual 0 -honorTolR 0
```

**Step policy (the driver's constants):**
- ds0 = 2e-5 m; × 2 after 6 good steps, up to 1e-3 m; ÷ 2 on a failed ladder.
- FLOOR at ds < 2e-7 m.
- Ladder: Newton (25 iterations) → NewtonLineSearch (40) → KrylovNewton (60, tolerance × 10).
- Test: `NormUnbalance` at 1e-5 × the applied vertical load (0.0415 kN).
- For determinism, 1 MKL thread with `MKL_CBWR=COMPATIBLE`, or `system Pardiso -deterministic` on current builds (WP-132).

**Why:**
- It reaches 1.74× the settlement of ModifiedEuler (0.0508 vs 0.0292) and is 1.3–1.45× cheaper per unit s/B (§8.2).
- Its per-increment error is controlled at 0.6–2e-4·p′ at every replayed checkpoint (§5.3). ModifiedEuler is ~20× worse only at the onset of plasticity and at parity from s/B ≈ 0.01. **Not "20× more accurate".**
- It **refuses** a state it cannot integrate, where ModifiedEuler commits ρ_α up to 13 and a spurious +6 % (§8.4).

**Rejected alternatives:**
- TolR 1e-3: earlier wall, no saving (§8.5).
- TanType 1: step-size collapse, 13× cost (§8.6).
- `-maxSubsteps` above 2000: only tested inside E_C2.

**What it does NOT buy: a capacity.** The wall stays. A SANISAND q–s past s/B ≈ 0.036 on this deck is carried through NonPosH refusals, and it ends at 0.0508 with no peak. Read q with the Krylov ±1.5 kPa band (§8.2).

## 13. Follow-ups — each is an owner/TIMs decision

1. **SAS-ME follow-up WP.** Carry IntScheme 129 as the campaign integrator (guide default, counters in the TIMs deliverables). Look at its two cost items on this deck: the maxSubsteps refusals on the surface ring (209 in E_B), and the re-seat chatter (10.1 M re-seats in E_B, per #892).
2. **Refusal-aware line search.** `rung_fail.csv` shows the LineSearch rung failing on a material refusal 119 times in E_B (41 in E_A, 180 in E_D). A line search that treats a refusal as "backtrack" rather than "rung failed" would keep those attempts on the ladder. Algorithm-side, no constitutive change.
3. **R1 — bounded h, opt-in.** h = b0 / max(a, c_A·√(2/3)·m), applied only where b:n ≤ 0. It is a **CONSTITUTIVE change** to DM04, so it needs an owner/TIMs decision before any code.
4. **R1b — hysteretic re-seat, a separate flag.** S4's cost (~27 M substeps per step) is the re-seat chatter; a re-seat hysteresis addresses that independently of R1.
5. **R2 — B/4, B/8, B/16 after R1.** Re-measure the mesh dependence of the wall and of the bands (§9) once the singularity is bounded.
6. **R3 — Perzyna viscoplasticity inside SANISAND,** only if TIMs' matched-settlement tolerance (ADR-90 OQ2) requires it. A nonlocal ψ̄ or a crack band does **not** treat this onset and is not proposed. No ADR-90 Duvaut–Lions regularization is recommended.

## 14. Verified / not verified

**Verified:**
- The gravity patch and resultants on every arm: K0 patch 2.4e-12 (B/16: 1.1e-11), resultant 1.1e-15.
- σ_zz recovery (1.5e-12).
- Self-replay of every dump by its own binary (≤ 3e-12).
- Determinism: A_r2 = A bit for bit.
- The local → Esmeralda handover: E_B = local B to 1e-5 kPa through step 40 (§4).
- The DP control reaches 0.15.
- The act's ModifiedEuler wall is reproduced (E_A, 0.0292).
- Every number in §0 and §8–§9 is recomputed from the committed `runs/E_*` records by `wall_table.py`, or taken from `esmeralda/analysis/tables/*`, which were computed from the same records plus the checkpoint fields on Esmeralda.

**Not verified here (quoted, with source):**
- The acoustic-tensor scan, the 10.1 M re-seats and the h-singularity mechanism: per #892.
- The 27× / +0.54 % TolR figures and the ±1.5 kPa Krylov band: campaign-verified, not recomputed on this branch.
- No material-point replay at the Esmeralda walls (§5.3).

**Not verified at all:**
- The ladders are interim (§11).
