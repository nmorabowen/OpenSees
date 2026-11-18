# WP-232k — soil models at MPM particles: method and results

ADR 232 §9 phase 5. Frozen interface: `tests/mpm_soil_interface.md` (D1–D23, with the round-2
amendments D19a and D21a). The kickoff facts, decisions K1–K9, gates K0 and S1–S10, the freeze,
and every ruling and owner decision are in `.claude/skills/ladruno-mpm/references/wp-briefs.md`
"232k". User guide: [[232_ladruno_mpm_guide]], section "Soil particles".

All gates were pre-registered before any soil-model MPM number. No band was widened. A miss is
recorded as a miss, with its diagnosis and the owner's decision. Units are kN, t, m, s and kPa;
stress is tension-positive.

## Method

- **What is admitted.** `LadrunoSANISAND` (IntScheme 1, 2 or 129; `-implex`) and `LadrunoHS`
  (and HSSmall), in 3D and plane strain. HS `-sig0`/`-sig03` and initial-stress offsets are
  refused (K3, Q2). Soil particles run serially (W14).
- **K0 restriction (owner).** SANISAND in mode F is admitted on B-spline bases only. `mpm analyze`
  refuses it on `gimp` and `linear`, and the message cites K0.
- **Staging.** One walk, `mpm stage <tag> <s>`. Each particle material's own `materialState` is
  set, so SANISAND's `Elastic2Plastic` runs at every particle. Three routes enter the walk:
  `mpm stage`, the `parameter … loadPattern` route, and mode A's `updateMaterialStage`.
- **Stage guard.** SANISAND's stage flag is class-wide, shared by every instance in the process. A
  guard refuses `mpm analyze` if that flag changed by a route that did not reach the particles.
- **Drivers.**
  - K0 and S2 use the column of D19.
  - S1 drives an affine grid motion (`mpm gridvel … affine`/`setAffine`): every particle sees the
    one-point FE record's dF.
  - S3 is the 232j L4 footing under K9's solution ladder.
- **Harnesses.** `tests/ladruno_mpm_soil_harness.py` (K0) and `tests/ladruno_mpm_soil_s_harness.py`
  (S1–S4, S10). Tests: `tests/test_ladruno_mpm_soil_*.py`.
- **References.** All are 232k-V's fork FE records under `tests/data/mpm/soil/` (never MPM).

## K0 — the SANISAND cell-crossing meter

The meter is m, the RMS over particles of |σ_p − σ_ref(Y_p)| / (ρgl0) over the xx, yy, xy and zz
components, after a 2.5 m (5-cell) rigid translation of the stressed column. The gate is
A = m(SANISAND) / m(DP) ≤ 2, or m(SANISAND) ≤ 0.01.

**The floor (amendment D19a).** The frozen floor was a gridvel box directly under the column. It
was a defect of the measurement:
- Each time a node row left the moving set, the soil's bottom row was left on GIMP tail weights.
  The full reaction was released at that row.
- Mode F elastic gave 0.188 against a twin of 0.0147. The bottom-row σ_yy error was −0.64 ρgl0
  and the peak inertia ratio 0.208.
- Mode A lost its support entirely.

The fix is a bonded elastic block under the column (excluded from the meter):
- Mode A prescribes the block's rows 2–4 cells below the soil.
- Mode F floats the block on an upward body force and translates everything by a uniform body
  acceleration. Local damping is off during the translation (D21a): FLAC damping is not Galilean
  invariant and would brake the rigid motion.
- The meter does not change with the settle length (0.5 s and 2 s agree to every digit).

**Mode A** (build `b94c9d75d`; `k0_meter_A_{gimp,bspline2}.json`; meter / twin):

| leg | gimp | bspline2 |
|---|---|---|
| elastic | 0.01408 / 0.01472 | 0.00232 / 0.00572 |
| DP (gated) | 0.01441 / 0.01505 | 0.00258 / 0.00663 |
| SANISAND | 0.01429 / 0.01472 | 0.00273 / 0.00572 |
| HS | 0.01356 / 0.01399 | 0.01734 / 0.02180 |
| dp_mildcone | 0.01812 / 0.01884 | 0.01790 / 0.02466 |
| **A** | **0.992 pass** | **1.058 pass** |

There were no post-flip floor events in mode A. SANISAND's drift was Δe ≤ 1.1e-4 and ‖Δα‖ = 0.

Mode-F records (esmeralda) are committed as `k0_meter_F_<basis>[_capped]_<leg>.json`; the
diagnosis probes as `k0_diagnosis/k0_F_gimp_capped_sanisand_<probe>.json`.

**Mode F, uncapped** (`b94c9d75d`, esmeralda). D23 was triggered: every translated soil leg had
post-flip floor events and no twin had any. These records are reported only:
- SANISAND: 0.1513 / 0.0148 (gimp), 0.0205 / 0.0055 (bspline2);
- HS: 0.0191 / 0.0142, 0.0173 / 0.0222;
- DP: 0.0179 / 0.0152, 0.0070 / 0.0066.

**Mode F, capped** (D23: a bonded one-cell elastic cap of 20 kPa; `15b109ad8`, esmeralda):

| | SANISAND | DP | A | verdict |
|---|---|---|---|---|
| gimp | 0.18557 (13 post-flip floor events) | 0.02335 | **7.95** | **miss** |
| bspline2 | 0.01616 (no floor events) | 0.02434 | **0.664** | pass |

Also capped: elastic 0.0335 / 0.0285, HS 0.0182 / 0.0206 (A_HS 0.78 / 0.85), dp_mildcone
0.0277 / 0.0311.

**Diagnosis of the gimp miss.** It is a numerical interaction, not a code defect. Records are
diagnostic only, on harness `dfb117de5`, and reproduce 0.18557 exactly.
- At the end of the post-flip hold the column is clean.
- The undamped explicit translation injects ±5–15 kPa of GIMP crossing noise.
- From t ≈ 9 s the two particle rows that share a cell split apart:

  | rows | p after the hold (kPa) | p after the translation (kPa) |
  |---|---|---|
  | Y 5.625 / 5.875 | 66 / 66 | 27 / 104 |
  | Y 7.125 / 7.375 | 47 / 47 | 2 / 96 |

- This is GIMP's grid-invisible, alternating two-row stress mode at ppc 2, and SANISAND locks
  it in. In that record alpha moved in 70 of 160 particles (max 0.58) and e by up to −0.0041, so
  plastic flow was engaged; whether the plasticity or a p-dependent elastic instability came first
  is not separated (critic C2).
- Probes:
  - Halving dt gives 0.112, so it is not a time-step effect.
  - Holding SANISAND in stage 0 gives 0.0335, the elastic value, so the plastic stage is the
    amplifier.
  - bspline2 gives 0.0162: the noise is ±3–5 kPa and the rows do not split.
- Engine paths checked: the particle stress equals the material's own stress exactly, and the
  stage walk's array refresh runs once, before the translation.

**Owner decision.** SANISAND is admitted in mode A on every basis, and in mode F on B-spline bases
only. HS is admitted everywhere. "Mode A on every basis" rests on quasi-static runs only.

**What the admitted passes measured (critic C1, verified).** In both mode-A records and in the
passing capped mode-F bspline2 record, SANISAND's back-stress never moved across the translation
(max |Δα| = 0 over 160 particles): the K0 stage-0 state sits on the cone axis and an oedometric
perturbation keeps it there. Those passes therefore measure SANISAND's pressure-dependent
hypoelastic response under cell crossing, not its plasticity. **K0p** (pre-registered in the AMEND
round) adds a stage-1 pre-shear (ν 0.05 and gravity +30 %) so that α moves before the
translation; its validity check (≥ 50 % of the particles moved α) runs before any meter. Windows
probe of the pre-shear on `7c932089c`: α moved in 160 of 160 particles in mode A gimp (max 0.085)
and in mode F bspline2 (max 0.159). The K0p meters wait for 232k-V's FE reference and run on
esmeralda, with K0 capped mode F on bspline3.

Mode-A A_HS = m(HS) / m(DP): 0.94 (gimp), 6.71 (bspline2), reported (HS's bspline2 twin is
0.0218 itself).

**Why the bspline2 twins carry a larger meter than the translated legs.** Mode A elastic on
bspline2 gives 0.0057 for the twin and 0.0023 translated. The twin's error sits at the soil–block
interface rows (Y 0.375: −0.028 ρgl0; Y 0.875: +0.009). At the end of the run the translated
column has the same alignment with the grid, but its rows show −0.003 / −0.001 there. The
B-spline quadrature error depends on where a particle sits within its cell:
- A particle that never moves accumulates the static error of one position.
- A particle that crosses five cells accumulates strain at many positions, so part of the
  position-dependent error averages out.

On GIMP the two are close (0.0147 against 0.0141), which fits GIMP's smaller position dependence.

## S rows

Builds `7b00fcf13` / `7c932089c`, Windows MSVC, unless noted.

| Row | Result | Verdict |
|---|---|---|
| S1 one-point identity (`s1_mpm.json`) | HS: 12 legs × 2 modes, max rel 1.2e-11. SANISAND: 22 of 24 pass (mode F on bspline2). SANISAND 3D mode A misses at step 12, the second plastic step, one step before the FE's own Gauss-point spread passes 1e-8: triax 3.2e-6, oedo 1.06e-8, with the FE's substep counts (5, 4). Option legs (`s1_onepoint_options.json`, rulings at f6ae63d79; build `f235a001f`): see "S1 option legs" below. | 22 pass; **2 pre-registered misses** (strict xfail, ruling at 926d27153); option legs: pass, except **IntScheme 2 shear, a pre-registered miss** (strict xfail, ruling d683e0db8) |
| S2 column + surcharge (`s2_mpm_hs_diagnosis.json`, `s2_mpm_sanisand.json`) | End of stage 0 and after the flip pass: HS 0.0128 / 0.0069, SANISAND 0.0127 / 0.0069. After the surcharge (registered top-row body force), max syy is HS 0.274 / 0.576 and SANISAND 0.265 / 0.776, at the top cell. | **recorded miss** (#1110; strict xfail) |
| S3 footing, HS (`s3_mpm_hs_r1_fbar.json`) | r1 fbar: q / q_FE 0.977 / 0.973 / 0.960 at s/B 0.01 / 0.02 / 0.05; reached the target; 0 failed attempts; 3.1 s. r1 off and r2 on reached on esmeralda; r3 and r2-off are running. | r1 pass; r2/r3 pending on esmeralda |
| S3 footing, SANISAND | Waits for the FE deck (probe v0 running; v1–v4 failed to reach). | pending |
| S4 stage parity (`s4_mpm_same_state.json`) | (i) The low-g column's flip map matches the rule exactly in both modes, with ≥ 1 reset particle each (g 0.25). (ii) Same-state leg: the **fallback branch** ran; see below. | pass |
| S5 refusal | SANISAND `-implex -maxSubsteps 1` aborts in both modes and in classic Tcl. Particles are unmoved bit for bit, the message names the material and tag, and later steps are refused until `mpm reset`. HS at γ = 10 fails loudly with no abort. | pass |
| S6 restart, mode A | SANISAND and HS mid-stage-1: the continued run is bit for bit. The restored host re-asserts SANISAND's stage flag (restart pack v4). | pass |
| S7 refusals and options | Symmetric-storage solvers refused (K7); unadmitted schemes refused; HS `-sig0` refused; Q2 offset refused; K2 (`ladrunoThreads 4` prints the serial warning, answer bit-identical); the K0 basis rule; `mpm info materials`. | pass |
| S8 unchanged answers | MPM suite: 615 passed, 16 skipped, 6 xfailed (all pre-registered), 0 failed; the same counts again on `f235a001f`. GCC on esmeralda: L4 r1 is bit for bit against the 232j record (346 of 346 values). | pass |
| S9 mutation | MPM family ZERO / SCALE: run by the orchestrator. Zone-A soil oracle: `test_s9_onepoint_triax` (HS and SANISAND, plane strain, both modes). | oracle in place; scores pending |
| S10 engines | Classic Tcl equals Python on the S1 hs_ps_triax leg (≤ 1e-12). | pass |

**S1, the 3D mode-A misses.** These are not kinematic errors:
- At step 12, mode A's particle strain differs from mode F's by 1.3e-12 relative. That comes
  through the SP displacement path; in mode F the particles agree among themselves to 7e-17.
- ModifiedEuler's adaptive substep sizes amplify the 1.3e-12 by about 1e6 while the substep count
  stays the same.
- Mode F holds 3.4e-10 at that step.
- The 1e-8 band presumes strain input free of round-off, which the FE's single brick has.

**S1 option legs** (`test_ladruno_mpm_soil_options.py`, build `f235a001f`, both modes, SANISAND
mode F on bspline2):
- **Pass under the S1 rule:** SANISAND IntScheme 129 and IntScheme 2 triax/oedo; HS `-implex`,
  `-smallStrain`, the HS deck's tension leg, and every `-sigmat` leg, the gated `tension` path
  included.
- **IntScheme 2 shear, a pre-registered miss (ruling d683e0db8).** Before the switch step (ps 21, 3D 61) every leg
  holds ≤ 1e-11. After it, the gate is the stress inside the FE's Gauss-point envelope widened by
  1e-3 relative. The record carries GP 1's stress and the GP spread, not per-GP min/max, so the
  envelope used is GP 1 ± spread, which contains [min, max]. ps mode A stays inside (worst 0.94 of
  the band). ps mode F has 26 of 110 post-switch steps outside (worst 1.46×, step 41); 3D mode A
  26 of 70 (1.71×, step 102); 3D mode F 9 of 70 (1.33×, step 102). Diagnosis: from the switch step
  the FE's 8 Gauss points take different CPPM branches (their substep counts part there). At step
  21 the ps particle sits at the far edge of the spread (distance to GP 1 = the spread to 1e-6).
  The particle's round-off differs from every GP's, so it follows a branch of its own; 8 GPs
  sample the branch set and do not bound it. Strict xfail; the miss goes to the owner at
  ready-flip.
- **SANISAND `-implex`: pass, compared at the material's committed response (ruling
  d683e0db8).** 232k-V's `-implex` records equal the scheme-1 records bit for bit. So does the
  MPM run at the material: the material's committed `stress` response and its 26-entry state equal
  the scheme-1 run's particle stress and state bit for bit, every step, all 12 legs and both modes
  (`test_option_implex_material_equals_scheme1`, the gated leg). The particle's own stress
  (`sigHat`) departs at step 11, the first plastic step, by 0.035 (ps oedo) to 0.233 (3D shear)
  relative. The engine keeps the stress the material returned at the trial call, which under
  `-implex` is the extrapolated stress the step was equilibrated with; `commitState` corrects the
  material, the particle keeps the trial value. Consistent with IMPL-EX, no C++ change; recorded
  as quirk `WP-232k-implex-particle-stress-is-trial` and a guide note.
- **`-implex` boundary-value leg** (`k0p_meter_A_gimp_implex.json`): K0p column, mode A gimp,
  capped, staging only, against `k0p_column_capped_implex.json` (232k-V, 185aacefd). α moved in
  160 of 160 particles (max 0.081); σ RMS 0.01862 ρgl0 against the gate 0.02. **Pass.**

**S1, rotation.** The FE record stores GP 1's material-frame stress, so S1 compares the particle's
unrotated σ̂ (`sigHat`). This is a reading of the record's convention, not a band change (ruling).
The record carries no FE Cauchy stress. The push-forward R σ̂ Rᵀ is therefore covered by the
slice-1 rigid-rotation objectivity case (ADR §8, "rigid-body translation and rotation of a
stressed body") and by 232j's L1(b) (a 30° rotation gives Δε ≤ 1e-12).

**S4 (ii), the same-state leg.** It took the fallback branch: the pre-flip states differ by
round-off, because mode A builds ΔF from displacements and mode F from velocities × dt.
- Row 0 (p < Pmin, the reset): pre-flip difference 2.2e-16, post-flip difference 0. Both modes
  commit σ = −Pmin I exactly.
- Row 39: pre-flip difference 3.4e-12, post-flip difference 3.4e-12, so the flip does not widen
  the gap.
- Both modes agree with the FE record to 6.7e-12 in σ and 1.6e-15 in state.

## Costs

| Run | Wall |
|---|---|
| K0 mode A | about 1 min per leg |
| K0 mode F | 5–55 s per leg on Windows; about 35 s on esmeralda |
| S1 | 0.1–7 s per leg |
| HS S3 r1 | 3 s |

SANISAND's M0b upper bounds (0.6–1.2 ms per plastic call) were not reached: increments at
particles are small.
