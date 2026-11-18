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
| bspline3 (`f235a001f`) | 0.00899 (no floor events) | 0.00543 | **1.656** | pass |

Also capped: elastic 0.0335 / 0.0285 / 0.00829, HS 0.0182 / 0.0206 / 0.02285 (A_HS 0.78 / 0.85
/ 4.21, reported), dp_mildcone 0.0277 / 0.0311 / 0.03809 (gimp / bspline2 / bspline3). On
bspline3 SANISAND's α did not move (‖Δα‖ = 0), as on the other capped passes.

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
and in mode F bspline2 (max 0.159).

**K0p results** (records `k0_meter_<mode>_<basis>_k0p_<leg>.json`, esmeralda, build `f235a001f`,
against 232k-V's `k0p_column_capped.json` at `after_hold2`, reference build `7c932089c`). The gate
is A_p = m(SANISAND) / m(DP) ≤ 2, or m ≤ 0.01; every leg is valid (α moved in 160 of 160
particles in the pre-shear) and has no post-flip floor event.

| mode, basis | SANISAND | DP | A_p | pre-shear max ‖Δα‖ | ‖Δα‖ in the translation | verdict |
|---|---|---|---|---|---|---|
| A gimp | 0.01822 | 0.01922 | **0.948** | 0.085 | 0.052 | pass |
| A bspline2 | 0.00607 | 0.03365 | **0.180** | 0.089 | 0.025 | pass |
| F bspline2 | 0.04727 | 0.04722 | **1.001** | 0.159 | 0.498 | pass |
| F bspline3 | 0.03215 | 0.04786 | **0.672** | 0.159 | 0.390 | pass |

Reported legs (elastic / HS / dp_mildcone): A gimp 0.0190 / 0.01747 / 0.01801; A bspline2
0.00366 / 0.0274 / 0.07715; F bspline2 0.01441 / 0.02734 / 0.10364; F bspline3 0.02867 / 0.03297 /
0.10619. In mode F the untranslated twins meter about as high as the translated legs (SANISAND
0.04634 / 0.04653 against 0.04727 / 0.03215), so the mode-F K0p meter is set mostly by the
pre-shear's dynamic state, not by cell crossing. α moves during the translation in mode F (0.50,
0.39) and much less in mode A (0.052, 0.025): the mode-F legs do measure SANISAND plasticity
under crossing.

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

Builds `7b00fcf13` / `7c932089c`, Windows MSVC, unless noted. The whole soil test set
(`tests/test_ladruno_mpm_soil_*.py`, `--runslow`), one pass on `f235a001f` with every record
committed: **160 passed, 10 xfailed, 0 skipped, 0 failed** (the xfails: S2 ×4, K0 capped gimp,
S1 SANISAND 3D mode A ×2, IntScheme 2 shear ×3).

| Row | Result | Verdict |
|---|---|---|
| S1 one-point identity (`s1_mpm.json`) | HS: 12 legs × 2 modes, max rel 1.2e-11. SANISAND: 22 of 24 pass (mode F on bspline2). SANISAND 3D mode A misses at step 12, the second plastic step, one step before the FE's own Gauss-point spread passes 1e-8: triax 3.2e-6, oedo 1.06e-8, with the FE's substep counts (5, 4). Option legs (`s1_onepoint_options.json`, rulings at f6ae63d79; build `f235a001f`): see "S1 option legs" below. | 22 pass; **2 recorded misses, owner-accepted 2026-11-18** (strict xfail, ruling at 926d27153); option legs: pass, except **IntScheme 2 shear, a recorded miss, owner-accepted 2026-11-18** (strict xfail, ruling d683e0db8; ps mode A on Linux only, owner-accepted the same day) |
| S2 column + surcharge (`s2_mpm_hs_diagnosis.json`, `s2_mpm_sanisand.json`) | End of stage 0 and after the flip pass: HS 0.0128 / 0.0069, SANISAND 0.0127 / 0.0069. After the surcharge (registered top-row body force), max syy is HS 0.274 / 0.576 and SANISAND 0.265 / 0.776, at the top cell. | **recorded miss** (#1110; strict xfail) |
| S3 footing, HS (`s3_mpm_hs_r<r>_<fbar\|off>.json`) | q / q_FE at s/B 0.01 / 0.02 / 0.05. Gated: r1 fbar 0.977 / 0.973 / 0.960 (0 failed attempts, 3 s); r2 fbar 0.969 / 0.968 / 0.907 (11 failed, 36 s, esmeralda), thin at s/B 0.05. Reported: r1 off 1.005 / 1.013 / 1.020; r2 off 0.987 / 0.992 / 0.982 (80 failed, 7 relaxed, 18 halvings); **r3 fbar 0.961 / 0.931 / 0.162** (56 failed, 10 relaxed, 5 halvings, 744 s), an artifact of the deck's footing condition (#1112), see below. | r1, r2 pass; r3 and off reported |
| S3 footing, SANISAND | No FE variant reached s = 0.1 m: v1–v4 (a15fdeb60) and v0 `tims_r1sep` (n4 stopped at the 2 h wall cap at s = 0.0215, q_max 636 kPa, 1341 relaxed steps; n8 at s = 0.0037). Per the freeze, reported only. Owner (2026-11-18): ONE MPM r1 leg with the v0 deck, fbar, gimp, mode A, K9 ladder, 2 h wall cap, reported beside the stalled FE (`s3_mpm_sanisand_r1_fbar_v0_tims_r1sep.json`, esmeralda, `f235a001f`): stopped at the halving floor at s = 0.0220 m (FE v0 n4: 0.0215 m); q(0.01 B) 586.3 kPa (FE 594.8, 0.986); 220 failed attempts, 41 relaxed, 46 halvings, 399 s; 31 SAS-ME refusals, all at the 2000-substep cap. See "S3 SANISAND" below. | reported, no reference |
| S4 stage parity (`s4_mpm_same_state.json`) | (i) The low-g column's flip map matches the rule exactly in both modes, with ≥ 1 reset particle each (g 0.25). (ii) Same-state leg: the **fallback branch** ran; see below. | pass |
| S5 refusal | SANISAND `-implex -maxSubsteps 1` aborts in both modes and in classic Tcl. Particles are unmoved bit for bit, the message names the material and tag, and later steps are refused until `mpm reset`. HS at γ = 10 fails loudly with no abort. | pass |
| S6 restart, mode A | SANISAND and HS mid-stage-1: the continued run is bit for bit. The restored host re-asserts SANISAND's stage flag (restart pack v4). | pass |
| S7 refusals and options | Symmetric-storage solvers refused (K7); unadmitted schemes refused; HS `-sig0` refused; Q2 offset refused; K2 (`ladrunoThreads 4` prints the serial warning, answer bit-identical); the K0 basis rule; `mpm info materials`. | pass |
| S8 unchanged answers | MPM suite: 615 passed, 16 skipped, 6 xfailed (all pre-registered), 0 failed; the same counts again on `f235a001f`. GCC on esmeralda: L4 r1 is bit for bit against the 232j record (346 of 346 values). | pass |
| S9 mutation | MPM mutation gate on CI (run 44318604262 at 42bbefac6): ZERO 0.500 (floor 0.46), SCALE 0.398 (floor 0.38). Zone-A soil oracle `test_s9_onepoint_triax` (HS and SANISAND, plane strain, both modes), demonstrated red by hand (`tests/mpm_soil_s9_red_probe.py`, record `s9_red_probe.json`, build `f235a001f`): unmutated, all 4 pass; skipping the stage flip goes red at step 11 (HS err_sig 0.036, err_state 0.999; SANISAND 0.014 / 0.41); skipping the material update on step 11 goes red at step 11 (HS 0.088 / 0.052; SANISAND 0.225 / 0.179); dropping the update's stress output goes red at step 2 (0.5), all against the 1e-8 band, in both modes. | pass |
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
  included. **HSSmall's cyclic shear** (the brief's R row) is covered here:
  `hs_smallstrain_{3d,ps}_shear` are `-smallStrain 140000 1.5e-4` on a 130-step, three-cycle
  shear path (cases block `s1_ref_onepoint_options`, `mpm_soil_cases.txt:266`; the FE's
  `-smallStrain` answer differs from the plain deck by 16 % at the end), and both pass the S1 rule
  in both modes.
- **IntScheme 2 shear, a recorded miss, owner-accepted 2026-11-18 (ruling d683e0db8).** Before the switch step (ps 21, 3D 61) every leg
  holds ≤ 1e-11. After it, the gate is the stress inside the FE's Gauss-point envelope widened by
  1e-3 relative. The record carries GP 1's stress and the GP spread, not per-GP min/max, so the
  envelope used is GP 1 ± spread, which contains [min, max]. ps mode A stays inside on Windows (worst 0.94 of
  the band); on Linux (Zone-A at ffec57a09) it falls outside from step 22, particle 11, the same
  branch mechanism on another libm/MKL round-off, so it is a strict xfail on Linux only:
  **owner-accepted 2026-11-18** under this ruling (merge on green; recorded on #1108). ps mode F has 26 of 110 post-switch steps outside (worst 1.46×, step 41); 3D mode A
  26 of 70 (1.71×, step 102); 3D mode F 9 of 70 (1.33×, step 102). Diagnosis: from the switch step
  the FE's 8 Gauss points take different CPPM branches (their substep counts part there). At step
  21 the ps particle sits at the far edge of the spread (distance to GP 1 = the spread to 1e-6).
  The particle's round-off differs from every GP's, so it follows a branch of its own; 8 GPs
  sample the branch set and do not bound it. Strict xfail, accepted by the owner (2026-11-18).
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

**S3, the r3 collapse: an artifact of the footing condition in the S3 deck (#1112; 232j's L4 r3
stop, #1107, sits at s/h = 0.75 too).** At s = 0.096 m
(s/h = 0.768, h = 0.125) q falls from 868 to 38 kPa in one step, then climbs like a fresh contact
(+54, +63 kPa per step); the run "reaches" s = 0.1 with q 155. It is neither physical (HS has no
softening; the FE n16 curve is monotone to 958 kPa) nor the ladder (a relaxed step still meets
NormUnbalance 1e-5, so it is equilibrated; the relaxed steps are a symptom). The mechanism:
- The footing is a set of prescribed grid nodes on the fixed row y = 0. `track y` only lowers the
  box's lower bound to y = −s, so the next node row joins at s = h.
- The soil under it moves down by about s. At ppc 2 the particle centres start at −h/4 and −3h/4,
  so the second row leaves the top cell at s ≈ h/4 and the top row at s ≈ 3h/4. Past that the cell
  next to the footing nodes is empty and the footing loses the soil until s = h.
- Both events sit at the same s/h on every rung. The first irregular q increment is at s/h = 0.256
  in r2 fbar (s = 0.064) and in r3 fbar (s = 0.032), with the same size (+10.5 against the usual
  +17.5); the collapse is at s/h = 0.768. r1 (s/h ≤ 0.2) meets neither; the off runs show no
  increment jump at 0.25 (F-bar's cell-centre average is what changes when a row leaves a cell).
- Trace (`tests/mpm_soil_s3_trace.py`, r2 run on to s = 0.2 m, Windows, `f235a001f`, records in
  `s3_diagnosis/`). fbar reproduces the esmeralda q (790.85 at s = 0.088). The footing's
  particles still in the top cell row drop 16 → 10 → 9 → 8 at s/h 0.248–0.264, on the first
  irregular increment. At s/h 0.760 the top row starts to leave (8 → 5), q drops from 1690 to
  1435 kPa, and then swings (1501, 1402, 1498, 1140) as top-row particles move back up and in and
  out of the cell. The off run, on Windows, stops at the halving floor at s/h 0.256, the same
  crossing (esmeralda's r2 off got past it with 80 failed attempts).

So s/B = 0.05 at r3 (s = 0.8h) is past the range in which this footing condition holds (about
0.75h); r1 (0.2h) and r2 (0.4h) are inside it. A footing that follows the soil (a node row inside
the tracked box before the top row leaves, or s_end/h ≤ 0.5 at every rung) is filed as #1112;
r3 stays a reported row.

**S3 SANISAND, the v0 deck (reported).** MPM r1 (fbar, gimp, mode A, K9 ladder) against the FE
v0 n4 record (`s3_fe_sanisand_v0_tims_r1sep_n4.json`, linearly interpolated at the MPM's s) over
the common range:

| s (m) | 0.002 | 0.004 | 0.006 | 0.008 | 0.010 | 0.012 | 0.014 | 0.016 | 0.018 | 0.020 |
|---|---|---|---|---|---|---|---|---|---|---|
| q MPM (kPa) | 96.8 | 155.2 | 210.8 | 265.2 | 319.1 | 372.5 | 426.0 | 479.4 | 532.9 | 586.3 |
| q FE (kPa) | 98.0 | 157.3 | 213.5 | 268.9 | 323.5 | 377.9 | 432.1 | 486.4 | 540.6 | 594.8 |
| MPM / FE | 0.988 | 0.987 | 0.987 | 0.986 | 0.986 | 0.986 | 0.986 | 0.986 | 0.986 | 0.986 |

The MPM accepted one more push (s = 0.022, q 639.9) before its halving floor; the FE stopped at
its 2 h wall cap at 0.0215 m. Both stop at the same depth on the same deck: the MPM shares the
FE's integrator stall (SAS-ME at its 2000-substep cap: 31 refusals, all `refCap`; 3.2e7 of 4.1e8
substeps rejected on error; 41 steps accepted on the relaxed rung, against the FE's 1341 of 1342),
so neither reaches s/B = 0.02.

**Granular-column collapse (R row, reported, not gated).** `tests/mpm_soil_collapse.py`: Lube's
column (aspect ratio 2, half width 1 m, symmetry on xmin, no-slip base), mode F on bspline2,
h 0.25 m, ppc 2 (128 particles), the k0_ref_column decks; dt 0.25 dt_c re-read every 0.05 s.
Recorded: runout (L∞ − L₀)/a from the front particle edge, H∞/2a, a history with max |v|, volume
checks (Σ V_p, det F) at the flip, the release and the end, floor hits after the release, flip
resets, the face q/p at the flip, SANISAND's Δα, and a refusal abort with its time.

*The first staging stood the soil columns up, for a reason in the materials.* Stage 0 elastic,
damped, on the unconfined column, then the flip and 4 s undamped (esmeralda, `f235a001f`): DP
fell (runout 3.41, H∞/2a 0.49) but SANISAND and HS ended where they started. The "runout −0.062,
height 0.969" in those records is the old metric (front particle CENTRE against a; 0.9375 and
1.9375 are the initial outer centres), not lost area: Σ V_p was 1.99990 (SANISAND) and 1.99978
(HS) of 2.0 at the flip and at the end, det F in [0.9997, 1.0013], and max |v| stayed below 0.016
m/s (SANISAND, front +0.0005 m in 4.5 s) and 0.0013 m/s (HS). Identical geometry because neither
moved, not shared state (each leg is a fresh process with its own material). The causes:
- SANISAND: the elastic column's free face carries q/p up to 2.8 against M_c 1.33, and vanilla
  ManzariDafalias's stage switch then RAISES M_c at every such point to contain the state ("stage-
  switch stress ratio exceeded the bounding surface; raising M_c", up to M_c 4.09 on a Windows
  rerun). A
  rate model with no driving strain keeps the stress; the column creeps.
- HS: it seeds its hardening from the committed stress at 0→1 and the deck carries c = 10 kPa: a
  vertical cut stands to 4c/γ · tan(45° + φ/2) = 4.2 m > 2 m. Standing is the right answer.
- DP: stage 2 returns the trial stress to its yield surface at the first plastic call, which
  unbalances the column.

*The gate staging (the harness default now).* Gravity and a smooth gate on the free face, both as
per-particle body forces (lateral pressure 0.5 ρg(2a − y) on the outer particle column); stage 0,
damped, 1 s; the flip under the gate (face q/p 1.83, no M_c raise); 0.5 s damped; release; 4 s
undamped. Windows, `f235a001f` (records `collapse/collapse_{dp,dp_mildcone,hs,sanisand}.json`;
they replace the first staging's; SANISAND's is the partial run the wall cap stopped, which a
6 h esmeralda rerun was not worth for a reported row):

| leg | runout | H∞/2a | Σ V_p at the end (of 2.0) | det F range at the end | floor hits | wall |
|---|---|---|---|---|---|---|
| DP (gated) | 3.49 | 0.52 | 2.454 | 1.000–4.77 | — | 13 s |
| dp_mildcone | 3.54 | 0.52 | 2.460 | 1.000–4.85 | — | 10 s |
| HS | 0.0008 (stands) | 0.999 | 1.9999 | 0.9994–1.0005 | 14 | 3 s |
| SANISAND | 0.81 at the 2 h cap (t = 3.25 s, 1.7 s after the release, still accelerating: max |v| 2.0 m/s) | 0.85 | 2.657 | 1.022–2.28 | 215 | 9011 s |

SANISAND now flows, slower to start than DP (front +0.08 / +0.20 m at 0.5 / 1.0 s after the
release, against DP's +0.61 / +2.11 m), then accelerating; Δα 0.010 under the gate and 8.3 after the release, the
large values at particles held at the p floor. Its cost in the flow is about 5 ms per particle-step
(128 particles, about 0.6 s per step), three orders above the K0p probe's: the flow's increments
are large. Two findings while building it: with M0b's Coulomb wall (μ 0.4) a base particle's centre
reaches y = 0 and the next step is refused as leaving the grid (A7; DP at t = 2.08 s), so the
base is no-slip; and a fixed dt set from the post-settle dt_c is refused by the CFL check during
the flow.

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

**Cost per particle-step at particles (R row; closes M0b's 222–920 range).** Probe
`tests/mpm_soil_cost_probe.py`, record `cost_probe_k0p_F_bspline2.json` (Windows, `f235a001f`):
the K0p column in mode F on bspline2, capped (160 soil + 32 floor-block and cap particles,
dt 5e-4), wall time
per 100-step chunk by phase, and every 20th chunk one extra step bracketed by F reads for the
strain increment a particle sees, |sym(F_{n+1} F_n⁻¹) − I|. The material cost is a leg's time
per particle-step minus the elastic leg's (same engine work).

| phase | median Δε per step (SANISAND; p90) | SANISAND | DP | HS | SANISAND / DP | HS / DP |
|---|---|---|---|---|---|---|
| pre-shear (α moves) | 1.7e-8 (1.5e-7) | 3.63 µs | 0.91 µs | 0.53 µs | **3.99** | 0.58 |
| translation | 1.5e-7 (5.6e-7) | 1.53 µs | 0.53 µs | 0.54 µs | **2.88** | 1.02 |

(per particle-step over all 192 particles; per soil particle × 1.2). An explicit particle sees
increments of 1e-8–1e-7 per step, three to four orders below M0b's 2.56e-4 and 9.44e-4, where
SANISAND / DP was 341 and 647 per plastic call. At these increments SANISAND costs 3–4× DP per
particle-step, not hundreds. The esmeralda K0p records agree in whole-leg wall time: (SANISAND −
elastic) / (DP − elastic) = 2.18 (F bspline2), 2.07 (F bspline3); mode A's legs run 0.4–3.2 s in
all and are too short to separate (9.7, 13.9).
