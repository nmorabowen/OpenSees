# WP-232k — soil particles: the interface (232k-C, round 1, the freeze candidate)

Status: **proposed** by 232k-C, 2026-11-17; binding once the orchestrator says "frozen". Authority:
the ADR (`Ladruno_implementation/232_ladruno_mpm_general_driver_adr.md`, §4–§6, §8, §9, amendment
"232k kickoff") and `wp-briefs.md` "232k" (K1–K9, K0, S1–S10). Anchors on `fc82c1a41`; lines drift.
232k-V's reference scripts were not read. Decisions are numbered (D*n*) for review.

## 1. K1 — admission

| | `LadrunoSANISAND` | `LadrunoHS` (and `LadrunoHSSmall`) |
|---|---|---|
| Type strings | `ThreeDimensional` → 33020, `PlaneStrain` → 33021 (`LadrunoSANISAND.cpp:2107`, `:2163`; the copy keeps the tag) | one class 33029, `dim_` from the string (`LadrunoHS.cpp:936-955`) |
| Hypo use | yes: the material gets the accumulated unrotated ε̂ (interface §9) | yes |
| Stage path | the K4 walk (§3) | the K4 walk; the lazy registry sync is complete on its own (§3) |
| Process statics | `mElastFlag` is **class-wide**, every SANISAND tag in the process, FE and MPM (`ManzariDafalias.h:557`); `wipe` zeroes it (`:560`) | the stage registry is **per tag**, FE and MPM instances of one tag share it (`LadrunoHS.cpp:54-73`) |
| WP-107 threads | no (`ManzariDafalias.cpp:1235`, delegated `LadrunoSANISAND.cpp:2758`) | no (inherits `Material::ladrunoThreadSafeUpdate`) |
| Length consumer | no (W9) | no |

- **D1 table.** `admitMaterial` (`LadrunoMPMPattern.cpp:391-395`) adds the three class tags, both
  modes. The existing order check (`:414`) holds: 6 / 3.
- **D2 SANISAND options.** Schemes 1, 2, 129 admitted, each with an S1 triaxial leg (mode A QS)
  against 232k-V's one-point driver. Every other scheme refused at admission ("IntScheme n is not
  admitted at particles: 1, 2, 129"); scheme 4 also carries function-scope statics
  (`ManzariDafalias.h:219`). `-implex` admitted with its S1 leg. Its default time source is
  `ops_Dt` (`LadrunoSANISAND.cpp:815`, "pseudo = ops_Dt (the default)"). Mode F sets it through
  `dom->applyLoad(t0 + dt)` (`LadrunoMPMPattern.cpp:767` → `Domain.cpp:2170` `ops_Dt = dT;`). Mode A
  sets it through the integrator's `applyLoad`, with Δλ = dt forced by the driver (implicit
  interface §12.1). A cut step changes `ops_Dt` exactly as it does in FE. `-implexControl`,
  `-meFallback`, `-flipAlphaIn`, `-implexFlipAbsorb` and the floors (`-Pmin`, `-Presidual`, `-pRe`, SAS-ME cutoff)
  are admitted unchanged (K5).
- **D3 HS options.** 3D and PlaneStrain only (other `dim_` refused by type). `-smallStrain`,
  `-sigmat` admitted. `-implex` admitted with an S1 leg: no time source, it extrapolates
  `dsig_p_n_` (`LadrunoHS.cpp:635-680`). `-sig0`/`-sig03` refused, §2.
- **D4 accessors** (fork files only, no vanilla row): `LadrunoSANISAND::ladrunoImplexEnabled()`,
  `static int ladrunoStageFlag()` (reads `mElastFlag`); `LadrunoHS::usesImplex()`,
  `hasInitialStress()`; `LadrunoSANISANDPlaneStrain::getStressZZ()` (= −`mSigma(2)`, as
  `getStressToRecord` `:175-182`). The material's own density is ignored at particles (mass is
  `-rho` of the body), as for DruckerPrager.

## 2. K2, K3

- **D5 K2: honoured as 1, with the warning that already exists. No new code.** Mode F has no
  threaded particle loop in v1, for any material ("the loop is serial (v1)",
  `LadrunoMPMPattern.cpp:77`). In mode A, `LadrunoMPMPoint` opts into nothing (`Element.h:143`, default
  `false`), so `Domain::update` runs the particles serially and prints once per domain "is not on
  the thread-safe allowlist … running the element loop SERIAL" (`Domain.cpp:2725-2731`).
  ElasticIsotropic particles get the same today. S7 tests that `ladrunoThreads 4` with a soil
  particle prints the warning and gives a bit-identical result.
- **D6 K3: HS `-sig0`/`-sig03` refused.** The stress is inside the material (yield sees it). But
  mode F's first step reads the host arrays: `fint = −Σ V σⁿ ∇S` from `P.sig`, which
  `addParticle` seeds from the Q2 offset (zero for soil), not from the material
  (`LadrunoMPMPattern.cpp:448`, core §6.3). The first step would run with σ = 0 while the material
  holds σ₀. The value is also uniform per tag, so it cannot carry a geostatic state. No gate needs it.
  Message: "LadrunoHS tag t uses -sig0/-sig03: refused at particles (K3: stage under gravity)".

## 3. K4 — one stage walk, both modes

Today: `MaterialStageParameter::setDomain` stops at the first answering element (`:79`,
`while (((theEle = theEles()) != 0) && (theResult == -1))`). Classic Tcl and Python both build it
(`TclUpdateMaterialStageCommand.cpp:61`, `OpenSeesNDMaterialCommands.cpp:367`). `LadrunoMPMPoint`
forwards only its own material (`LadrunoMPMPoint.cpp:1108-1115`), so in mode A one particle gets
`Elastic2Plastic()` (`ManzariDafalias.cpp:975`). The others flip lazily in `ladrunoTrialUpdate`
(`LadrunoSANISAND.cpp:3397-3406`), and that path runs `ladrunoRunStageFlipOnce` (`:5054-5136`)
without the p < 3·Pmin reset (`:5789-5794`).

A trap: `parameter … updateMaterialStage` + `updateParameter` arrives as a **double** in both
engines (`TclParameterCommands.cpp:327-334`; `OpenSeesParameterCommands.cpp:498-504`). SANISAND
id 1 reads `info.theInt` (`ManzariDafalias.cpp:973`), which is then indeterminate
(`LadrunoSANISAND.cpp:5253-5257`). HS id 1 has the same trap (`LadrunoHS.cpp:1227`).

- **D7 the walk.** `LadrunoMPMPattern::stageWalk(matTag, stage)` builds a host-local `Parameter`,
  never added to the Domain. Every particle material with that tag registers through its own
  `setParameter({"materialState", tag})`: SANISAND id 5, HS id 13, DP `materialState`. Then
  `update(double(stage))` runs. That is FE's directly reached Gauss-point path, so every SANISAND
  particle runs `Elastic2Plastic()` and then `ladrunoRunStageFlipOnce` (`:5309-5313`, which sets
  `mFlipSeen`). `mats[p]` is the material in both modes: mode F owns it, and in mode A it "stays a
  view" (`LadrunoMPMImplicit.cpp:521`).
- **D8 host arrays.** `Elastic2Plastic` writes `mSigma = mSigma_n = Pmin·I` at low p. The walk
  rewrites `P.sigHat` and `P.sig` (push with `P.R`) from `getStress()` for every particle whose
  material stress changed bitwise. Without that, mode F's next P2G would use the pre-reset σ.
  DruckerPrager answers stay bit for bit (S8).
- **D9 the routes, all into D7.** (a) `mpm stage <matTag> <stage>`, new, both modes, the
  documented route. (b) `parameter <p> loadPattern <mpm> materialState|updateMaterialStage <matTag>`:
  for these names the host registers **itself** instead of forwarding (`:347-358` today forwards).
  Its `updateParameter` walks with `(int)theDouble`, which also removes the theInt trap.
  (c) Mode A `updateMaterialStage`: `LadrunoMPMPoint::setParameter` sees
  `param.getClassTag() == PARAMETER_TAG_MaterialStageParameter` and hands the host the
  registration (once per Parameter). The host walks with `theInt`, which is the int path both
  commands use. The first-answer stop is harmless when the first answerer is a particle.
- **D10 the guard, so no particle can miss the reset.** At every `mpm analyze`, before any trial,
  the host compares `LadrunoSANISAND::ladrunoStageFlag()` with the stage it last walked (0 at
  admission and after `mpm reset`). Any mismatch refuses `mpm analyze`, naming the tag:
  "LadrunoSANISAND's class-wide stage is 1, the particles were staged 0 — `mpm stage <t> 1`". A
  mismatch comes from an FE element answering first, from mode F's `updateMaterialStage` reaching
  an FE element (W12), or from another SANISAND tag's flip. A particle's lazy path runs only in a
  trial update, and trials run only inside `mpm analyze`: mode-A elements refuse `update()`
  outside it (implicit interface D1.9). So the guard runs before the lazy path can.
- **D11 HS needs no guard.** Its lazy path is the full transition: `syncStageFromRegistry` →
  `applyStageTransition` (`LadrunoHS.cpp:503-513`), the function that `updateParameter` 1/13 runs
  (`:1225-1232`). It seeds from the committed `st_n_.sig` and keeps the stress (`:517-577`).
- **Warning (guide):** SANISAND's stage is class-wide. A SANISAND FE model and SANISAND particles
  in one process share it, and an FE mesh flipped by `updateMaterialStage` resets one element's
  points only. For 232k-V: FE column and footing points flip lazily without `Elastic2Plastic`. That
  is inert on the gated domains: K0 stress ratio η = 3(1−K0)/(1+2K0) = 0.86 < M_c = 1.33, and the
  top particle's p ≈ 1.6 kPa > 3·Pmin = 0.30 kPa. Only S4's low-p particles see the difference.

## 4. K6 — refusals (an objection: neither proposed trigger refuses at commit at particles)

Both proposed triggers are **trial** refusals, and MPM forwards trial codes. In mode F the step
fails and nothing aborts: "failed to integrate the trial strain", with `revertAll`
(`LadrunoMPMCore.h:990-995`). In mode A the element returns −4 (`LadrunoMPMPoint.cpp:339`), and
the step is cut and retried (implicit interface §3.8). HS `|Δε|² > 25` (`LadrunoHSKernel.h:1721`)
and the SAS-ME start refusal both take this path. SANISAND's commit latches
(`LadrunoSANISAND.cpp:4884-4926`) fire only under an element that **discards** the trial code. HS
refuses at commit only through `refused_` (`LadrunoHS.cpp:831-834`), which `revertToLastCommit`
clears (`:871`). So HS has no reachable W17 path at particles.

- **D12 SANISAND, the W17 abort:** `-implex -maxSubsteps 1` (1 is the smallest cap; 0 means
  uncapped, `:248`), then one stage-1 affine shear step Δγ = 1e-3 on the S1 block. The companion
  hits the cap in `commitState` and declares `ladrunoNoteCommitRefusal()` (`:4623-4683`). Mode F:
  `pendingRefusals` → −2 → `aborted`, particles unmoved (`LadrunoMPMCore.h:1008-1016`,
  `LadrunoMPMPattern.cpp:792-798`). Mode A: batched G2P skipped (`LadrunoMPMImplicit.cpp:135`) →
  −4 → `aborted` (`:1224`, `:1255`). Round 2 first checks on the one-point FE that the trigger
  fires.
- **D13 naming the first refuser:** the particle commit is serial in both modes. `Mats::commitState(p)`
  and `LadrunoMPMPoint::commitState` record the first p whose commit raised
  `ladrunoPendingCommitRefusals()`. The message adds "particle id, nDMaterial class and tag".
- **D14 HS (amend S5's HS leg):** `|Δε|² > 25` from one affine simple-shear step γ = 6
  (Hughes–Winget leaves γ unbounded, `LadrunoHypoKernel.h:276-330`), with `mpm cut 2 0` in mode A.
  The gate: the step fails loudly, the message names the material and tag, the particles keep
  the last accepted step, and there is **no abort**: a later `mpm analyze` is admitted.

## 5. K7, K8

- **D15 K7:** 232j's C3 list (`LadrunoMPMImplicit.cpp:406-425`) refuses when `antiLocking`, **or**
  when any particle material is SANISAND or HS without `-implex`. Under `-implex` both tangents are
  frozen isotropic elastic, so symmetric (`LadrunoHS.cpp:669`; SANISAND
  `ladrunoImplexFreezeTangent`). The message names the material.
- **D16 K8: no contract change.** New queries only. `mpm particleData sigZZ`: in plane strain,
  `getStressZZ()` of each particle material (ElasticIsotropic, DP and HS have it,
  `ElasticIsotropicPlaneStrain2D.cpp:138`, `DruckerPragerPlaneStrain.cpp:122`,
  `LadrunoHS.cpp:779`; SANISAND gets D4). In-plane R leaves σ_zz invariant. `mpm particleResponse
  <args…>`: the material response per particle, flat in storage order. `mpm particleResponse -size
  <args…>` gives its length. SANISAND `state` is 26 values: α at 6–11, e at 24
  (`ManzariDafalias.cpp:1072-1083`). HS has `voidRatio`, `p_p`, `gamma_p`. `mpm info materials`: per
  tag, the class, scheme, implex, stage, count (S7).

## 6. The K0 meter (two objections, each with a fix)

- **D17 `track` cannot carry a floor.** It moves only the lower bound
  (`LadrunoMPMImplicit.cpp:661`, `lo[gv.track] += gv.v[gv.track] * std::max(t - gv.start, 0.0);`). A
  floor under the column needs its top to descend, and before `start` the gridvel is absent, so
  the column would fall during staging. **Fix:** `mpm gridvel … profile <t0> <T_r> <t1> translate
  <axis>`, present from t = 0. The velocity is v·φ(t), with φ a trapezoid (0 before t0, a linear
  ramp over T_r, 1, a ramp down that ends at t1, 0 after). `translate` moves both bounds by v·Φ(tₙ)
  at the start of step n, with Φ = ∫φ. Mode A sets each step's SP value from
  v·(Φ(tₙ₊₁) − Φ(tₙ))/dt (`LadrunoMPMPattern.cpp:215-222`). `start`/`track` are unchanged.
- **D18 mode F has no gridvel at all.** `mpm gridvel` and `mpm pforce` are accepted and then
  ignored: the core has no prescribed velocity (`LadrunoMPMCore.h:391` `FaceConditions` only), and
  `analyze()` never reads `gridVels` (`LadrunoMPMPattern.cpp:690-800`). That is a silent wrong
  answer today. **Fix:** a prescribed velocity in step 4's full operator and in the
  second-projection operator. On the start-of-step node set, the component is
  v·(Φ(tₙ₊₁) − Φ(tₙ))/dt. A mode-F `pforce` is refused loudly. Test: an elastic block carried by a
  `profile` floor moves exactly v·Φ, and a no-gravity block gets zero stress.
- **D19 the case**, kPa / t / m / s (ρ = 2.0 t/m³ is the brief's 2000 kg/m³; g = 9.81;
  ρgl0 = 196.2 kPa). **Two cells wide**: a B-spline grid "needs at least 2 cells on every axis"
  (`OPS_LadrunoMPM.cpp:195`). That is an amendment to "one cell wide"; the problem stays 1D.
  Column x ∈ [0, 1.0], y ∈ [0, 10], ppc 2 (160 particles). `origin 0 -3.5`, `h 0.5`, `cells 2 28`.
  `fix xmin x`, `fix xmax x`, and **no `fix ymin`**. Floor:
  `gridvel floor box 0 -1.0 1.0 0.0 v free <−v> profile <t0> <T_r> <t1> translate y`. Its x range
  stays on real nodes, which is B10-safe on B-splines (`LadrunoMPMPattern.cpp:669-676`). The
  translation is D = 2.5 m (5 cells). At the end the floor box is [−3.5, −2.5] and the lowest
  particle is at −2.375. Neither support (≤ 1.5h) reaches the ymin ghost row at −4.0, so ghosts
  and `fix ymin` cannot interfere. The walls are slip faces with ghost layers (Amendment 5 F1,
  sliding along a face). Twin: v = 0, the same timeline.
- **D20 mode A:** QS `LoadControl`, `constraints Transformation`, which a non-homogeneous gridvel
  needs (`LadrunoMPMImplicit.cpp:385-388`), UmfPack (D15), NormDispIncr 1e-10. Gravity is
  `series linear 1.0` with dt 0.1 (10 steps), then `mpm stage`, then a 2-step hold. The
  translation runs at v = 1 with dt = 0.0625 = h/8 (40 steps, T_r = 0), then a 2-step hold: 54
  steps.
- **D21 mode F:** damping 0.8, flip 1, dt = 5e-4 s. SANISAND's stage-0 M = K + 4G/3 = 300 MPa with
  G = G0·P_atm·(2.97−e)²/(1+e) = 81.6 MPa. The stage-1 factor √(p/P_atm) ≤ 1.14 gives c ≤ 415 m/s
  and dt_c·s = 6.0e-4. Timeline: gravity ramp 1.0 s, settle 0.5 s, `mpm stage`, hold 0.25 s. Then
  translation with v = 0.1 m/s, T_r = 2 s, t1 = t0 + 27 s, then settle 0.5 s: 58,500 steps.
  **Inertia check:** the ramps give ρ·l0·v/T_r = 1.0 kPa = 0.51 % of ρgl0. The velocity is
  continuous, so no ρcv jump: a constant v = 0.1 jump would give 83 kPa = 42 %, and a jump within
  1 % would need v ≤ 2.4e-3 and 2.1e6 steps. Measured as |dP_y/dt| / W from `mpm info summary`,
  averaged over 100 steps: gated ≤ 1 % at the measurement, maximum reported. FLAC local damping is
  not Galilean invariant: under drift `sign(v)` stays fixed and it cannot dissipate
  (interface §6.5). So the meter is read after the ramp-down and settle, not during the drift.
- **D22 harness reads:** Y_p is `particleData x` at step 0. σ is `particleData sig` and `sigZZ`.
  Crossings come from x every 10 steps (one crossing takes ≥ 10⁴ steps in mode F). Δe and ‖Δα‖
  come from `particleResponse state` against the twin. Floor, clamp and refusal counts come from
  the material responses `commitStats` / `sasStats`.
- **D23 K5 contingency, pre-registered:** the K0 top is a free surface. The top particle has
  p ≈ 1.6 kPa, 5× above 3·Pmin. If a gated K0 run records any floor or clamp particle-step, that
  run is not admitted. The row then re-runs with a bonded elastic cap: one cell thick, 20 kPa, the
  FE reference with it. That rule is fixed now, before any number.

## 7. Costs (upper bounds from M0b: SANISAND 0.61 ms per plastic call at Δε 2.6e-4, 1.16 at 9.4e-4)

| Run | Material calls | SANISAND | DP / elastic | HS |
|---|---|---|---|---|
| K0 mode A (×16) | 54 steps × ~6 it × 160 = 5e4 | ≤ 1 min each | seconds | ≤ 1 min |
| K0 mode F (×16) | 58,500 × 160 = 9.4e6 | ≤ 1.6 h each (≈ 5 min if a tiny-Δε call is ~30 µs) | ≈ 30 s | ≤ 10 min |
| S2 (×4) | ~25 steps × 6 × 160 | < 1 min | — | < 1 min |
| S3 r1 to s = 0.05B | DP's r1 (88 failures, 24 s to 0.2 m): ≈ 3e6 | 0.5–1 h | — | ≈ 10 min |
| S3 r2 | DP's r2 (374 failures, 450 s): ≈ 4.6e7 | **8–15 h** | — | ≈ 1–2 h |
| S3 r3 | 4× r2 per iteration, DP stalled at s = 0.094 m | 30–60 h: **over 12 h, not run** | — | — |

On esmeralda the K0 mode-F runs go side by side: about 1.6 h of wall time. These are upper bounds:
the 232k R row measures the real increments.

## 8. Cases and tests

`tests/mpm_soil_cases.txt` is created by 232k-V (its reference blocks). This worker adds the MPM
blocks. New keys: `material <t> sanisand <deck…> [opts]` and `hs <…>`;
`stage <t_pseudo> <matTag> <s>`; `gridvel … profile … translate`; `gridvel <n> box … affine
<L ndm×ndm> [ref <x>]`. The affine key drives S1: every node is prescribed in mode F, the
perimeter in mode A on a 2×2(×2)-cell block with the interior free, so mode A keeps equations to
solve. The harness sets L per step from 232k-V's ΔF path (`mpm gridvel <n> setAffine <L>`). New
observables: `meter_rms`, `crossings`, `state_drift`, `inertia_ratio`, `floor_counts`. Test files
`tests/test_ladruno_mpm_soil_*.py`:

| File | Rows | Zone-A |
|---|---|---|
| `_admission` | K1, D6, D15, S7, K2 warning | yes, < 10 s |
| `_stage` | D7–D11, S4 on a 2×4-cell column, both modes, both engines | yes, < 20 s |
| `_refusal` | S5 (D12–D14), classic Tcl and Python | yes, < 15 s |
| `_onepoint` | S1 (S9's small triaxial and S10's engine parity in zone_a, the rest marked slow) | ≈ 30 s |
| `_gridvel` | D17–D18 profile, translate, mode-F floor | yes, < 10 s |
| `_restart` | S6 | yes, < 15 s |
| `_k0`, `_column`, `_footing` | K0, S2, S3: they re-check the JSON from the warrant run; the runners are esmeralda-only | JSON check < 2 s |

S1's FE identity: no FE element drives `hypoIncrementFromF`. `-geom hypo` (LadrunoBrick, 3D only)
drives `hypoIncrement`. The two have the same algebra but different rounding
(`LadrunoHypoKernel.h:254-275`). So 3D S1 is held to rel 1e-8 against a one-element LadrunoBrick
`-geom hypo`, and plane strain takes the small-strain rule (‖ε‖ ≤ 5e-4, rel 2e-3).
