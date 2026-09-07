---
title: ADR 94 review verdict — ASDPlasticMaterial3D
project: Ladruno
status: verdict shipped — fixes are separate WPs
owner: nmora
tags: [implementation, material, review, verdict]
---

# ADR 94 review verdict — `ASDPlasticMaterial3D`

Build `52314165a` (== HEAD of the review branch at R0; `ops.ladrunoBuild()` verified).
Plan: [[94_asdplastic_review_plan]] (merged #802). Execution PR: #804. Lane files:
`_adr94_inventory.md`, `_adr94_hlist_R1A/R1B/R1C.md`, `_adr94_hb_drift.md`,
`_adr94_redblue/{red1_cpp,red2_numerics,red3_tests_process,blue}.md`,
`_adr94_components.md`, `_adr94_matrix.md`, `_adr94_contract.md`. Tests:
`tests/test_adr94_*.py` (9 files). Owner decision D3 (2026-09-06): SP/MP out of scope.

## §0 Verdict

**Not warrant-grade as shipped; adequate for the workloads that have actually used it,
provided three conditions hold.** Every result the fork has published on this material
(Cerro Lindo M-series, the ADR-84 battery, the finite-strain D3 oracle) was obtained on
`Backward_Euler`, on Mohr–Coulomb-family yield functions with `H = 0`, in kPa, on
compression-dominated paths. Inside that envelope the converged stresses are exact
(blue's two-cube test, R1-A's VM closed form at 3.7e-14, R1-C's HB compression parity)
and the defects below cost iterations, not answers. Outside it — softening, shear-dominated
VonMises/Drucker–Prager, rock in net tension, any non-default integrator, any host other
than the tet — the material can commit inadmissible or NaN states while `analyze()` returns
0. Nothing in the shipped code can be made to fail loud on `stdBrick`, and only two of
fifteen failure sites can fail loud on `LadrunoBrick`.

Sixteen hypotheses were tested (the plan's H1–H15 plus the convention split R3a found).
Fifteen are CONFIRMED in whole or part, one (H6's accuracy half, on VonMises) is REFUTED,
and H3 was scoped out by D3. Twenty-two additional findings came from the red/blue and
component lanes. The fix list (§7) is split into fork-local opt-in work and jaabell-bound
framework work, per the ADR-84 §3 bit-identical strategy.

## §1 Findings register (adjudicated, severity-ranked)

Severity is for the fork's real use, after adjudicating red against blue (§2). "RESULTS"
means a converged step can commit a wrong or inadmissible state; "COST" means converged
states are exact and only Newton effort or diagnostics suffer.

| # | Finding | Sev. | Effect | Pinned by |
|---|---|---|---|---|
| **B1** | **H13** Unknown integration-option and model-parameter tokens are silently swallowed; unset parameters default to 0 (a typo'd `MC_phi` runs at φ = 0). | blocker | RESULTS, silent | `test_adr94_hlist_mechanical::test_H13_*` (2) |
| **B2** | **Fail-loud is unreachable on hex hosts** (red-numerics Q1, R5). `LADRUNO_MATERIAL_REFUSED` is returned at 2 of 15 failure sites, both in `Backward_Euler`; the other 13 return a bare −1 that `LadrunoBrick` drops (sentinel-only compare, all 6 paths) and `stdBrick` drops along with everything else. Same deck: `LadrunoBrick` 0/20, `stdBrick` 20/20 all zero. | blocker | RESULTS, silent | `test_adr94_redblue_numerics::test_R2_strict_convergence_is_a_noop_on_stdbrick`, `test_adr94_contract` (8) |
| **B3** | **H7** The `dLambda + deltaLambda < 0` branch commits the elastic predictor exactly, returns 0, and `strict_convergence` does not gate it (an eighth silent-accept site in the default integrator). Fires whenever softening exceeds `n:E:m` (`H_iso = −120000` vs `2G = 53846`): `f_VM` grows +3.5 → +87.5 kPa over four "converged" steps. | blocker (softening decks) | RESULTS | `test_adr94_hlist_numerics::test_H7_*` (2) |
| **B4** | **H10b + red-numerics Q5** Drucker–Prager hydrostatic tension through the apex commits **NaN with `analyze() == 0`**. Root cause is two-part: the NaN originates in `VoigtVector pressure_part; pressure_part *= 0.0;` on uninitialised Eigen storage (`DruckerPrager_YF.h:64-65`, `DruckerPrager_PF.h:68-69`, the ADR-84 constructor trap again; also `NullHardeningTensorPolicy` and AF's saturation branch), and the BE NaN guard *does* fire with a bare −1 that the host drops (B2). The apex-return body in BE (2093–2161) is commented out for every YF; `RoundedMohrCoulomb_YF` declares the apex trait without defining either method. | blocker (DP in tension) | RESULTS, NaN | `test_adr94_hlist_hb::test_H10_dp_apex_*`, `test_adr94_redblue_numerics::test_R2_dp_apex_nan_guard_fires_and_is_swallowed` |
| **B5** | **H16 (new, R3a)** Shear-slot derivative convention is split. `VonMises_YF`/`_PF` return the tensor derivative (shear slot = ∂f/∂σ12); MC, HB, MCTC, StiffSoil return the Voigt derivative (2×); `DruckerPrager_YF` matches neither (0.97 relative error on normal slots, a genuine gradient error). Consumption is split the same way: BE contracts with `tensor_dot_stress_like` (doubled shear), the other five sites with a plain Voigt dot, and every site accumulates `dLambda*m` and `Eelastic*m` expecting engineering-shear `m`. Net: VM flow direction under-counts shear at 5 of 6 sites and in the plastic-strain update; MC-family over-counts shear in BE's consistency scalar (convergence rate only, since BE iterates to Φ = 0). | major → blocker-candidate | RESULTS on shear-dominated VM/DP paths (magnitude **unmeasured** — the runtime probe did not converge in budget); COST for MC-family in BE | `_adr94_components.md` "Convention adjudication" (analytic, exact); `test_adr94_components` (VM Continuum-vs-numerical 1.4–1.8 %) |
| **M1** | **H1 + F1** Class-static `Stiffness` (and `dsigma`, `depsilon_elpl`, `intersection_*`): every GP, element and tag of one specialization is assembled with the tangent of the last GP integrated; `getInitialTangent()` mutates it as a side effect. Two disconnected elements with tangents 13.6 % apart assemble bit-identical blocks. Converged results are exact (blue: 1e-6 agreement with each cube solved alone); cost is 62 % more iterations on that model and non-convergence on harder ones. Permanent ADR-75b threading blocker; the YF/PF `static vv_out` buffers (F2) widen the sharing to every combo that reuses a functor type. | major (COST), blocker for threading | COST; RESULTS only via cutback history | `test_adr94_hlist_numerics::test_H1_*` (3), `test_adr94_redblue_blue` (2), `test_adr94_redblue_cpp` |
| **M2** | **H5 + red Q1** `strict_convergence` reaches only `Backward_Euler`; FE, FE_sub, ME, RK45_old commit `f_MC` in the hundreds–thousands with the flag on; BE_LS and RK45 unguarded by reading. Four more accept-without-check sites: ME/RK45 accept unconditionally at `dT_min` (and clamp `dT ≥ dT_min`, so a hard step always reaches it); `One_Step_Return` never verifies itself and uses a discarded stage's `Eelastic`. | major | RESULTS on non-default integrators | `test_adr94_hlist_mechanical::test_H5_*` (4 cases) |
| **M3** | **H6 (tangent half)** No `tangent_type` reproduces the consistent tangent of the committed map (`Continuum` 57 %, `Secant` — the default — 80 %, `Elastic` 103 %, `Numerical_*` 31 %, which differentiate a third map). Default costs 5.3× `Continuum`'s iterations. The map is Ortiz–Simo cutting-plane, not closest-point; its IV update is a Newton-path quadrature, exact only for constant `h` (VM + normalised flow); AF's recovery term is integrated explicitly inside the implicit return. **Accuracy half REFUTED on VM** (3.7e-14 vs closed form; normal does not rotate). | major | COST (VM, MC perfectly plastic); RESULTS for AF / rotating-flow hardening (latent, 22/46 specializations carry AF) | `test_adr94_hlist_numerics::test_H6_*` (3) |
| **M4** | **H10a** `HoekBrown_YF` was one commit behind `jaabell/ASDP` (`60d9b9b23`): our if/else was discontinuous at `arg = 0`; a uniaxial tension path locked onto `σci·s` = 587 kPa instead of `σt` = 245 kPa (factor `mb` = 2.4) then stalled. Compression paths bit-identical between trees. **shipped wp/94d** — jaabell's composite `max(f_shear, f_tension)` ported wholesale; tension now yields near the textbook 245 kPa with no stall; compression rows unchanged. | major (rock in tension) | RESULTS | `test_adr94_hlist_hb::test_H10_hb_*` (3) |
| **M5** | **Tolerance is absolute in stress units** (red Q2). Default `f_absolute_tol 1e-6` against `|Φ|` scaled by σy / c·cosφ / σci·sᵃ: the same MC problem passes 20/20 in kPa and is refused on step 1 in Pa. HB's central-difference normal through the discontinuity adds ~3e6 of spurious gradient inside `HB_ds`. | major | RESULTS (refusal vs acceptance decided by units) | `test_adr94_redblue_numerics::test_R2_f_absolute_tol_makes_strict_convergence_unit_dependent` |
| **M6** | **H4 + F4** `revertToLastCommit()` is a no-op; `revertToStart()` returns −1 that `Domain::revertToStart()` and `OPS_resetModel` ignore. `ops.reset()` leaves the material inconsistent with geometry (a third stress value, neither zero nor the pre-reset commit). Cutback-recovery divergence measured at ~6e-9, inside Newton tolerance (inconclusive). | major | RESULTS after `reset`; latent on fixed-step decks | `test_adr94_hlist_mechanical::test_H4_*`, `test_adr94_contract` (revert tests) |
| **M7** | **H8** `Backward_Euler_LineSearch`: hardcoded `max_iter = 30`, ignores `n_max_iterations` and `strict_convergence`, the "line search" tests a linear prediction so α = 1 always passes, the split loop solves one reduced increment and reports success for a strain the element never asked for. 2/20 vs BE's 20/20 on the ADR-84 leg. | major (if selected) | RESULTS if selected | `test_adr94_hlist_numerics::test_H8_*` (3) |
| **M8** | **H9** ME and RK45 drift checks are empty `if` blocks; ME commits `f_MC ≈ 1.1e5` with `return_to_yield_surface Disabled`. `Runge_Kutta_45_Error_Control_old` is still parse-reachable; `Full_Backward_Euler` and the Crisfield/Multistep enum values are accepted by the setter with no dispatch case. | major (explicit integrators) | RESULTS | `test_adr94_hlist_mechanical::test_H9_*` |
| **M9** | **H15** BE and BE_LS evaluate `E(σ)` once at commit; ME/RK45 per stage; `One_Step_Return` uses a discarded stage's `E`. StiffSoil/DuncanChang are different materials per integrator. A `StiffSoilShear` triaxial drive NaNs on step 1 for every parameter set tried; `StiffSoilShear_PF` is non-finite at 6/194 cloud points. | major (StiffSoil) | RESULTS for stress-dependent elasticity | `test_adr94_hlist_mechanical::test_H15_*` (structural), `_adr94_components.md` |
| **M10** | **F3** `sendSelf`/`recvSelf` print "not implemented" and return 0. Out of scope by D3; recorded because `database`/`save` paths are serial and also affected. | major (scoped out) | RESULTS under MP/database | none (D3) |
| **m1** | **H12 + AF** 161 `cout` vs 9 `opserr`; commit prints one line per material per step; AF prints 6–7 lines per IV per iteration per GP; per-GP `std::map::operator[]` option lookups. | minor | COST, diagnostics | `test_adr94_hlist_mechanical::test_H12_*` |
| **m2** | **F5** `setParameter` ids 7/8/16–21 write `CommitStress` only but set `stress_set_externally`, bypassing the `InitialP0` seed; the tag-match guard is `if (true)`. | minor | RESULTS for staged/geostatic decks using those ids | `test_adr94_redblue_cpp` |
| **m3** | **H2** `getClassType()` returns a dangling `c_str()` (SSO luck). **H14/F6** `getCopy()` omits `first_step` (latent: hosts copy at construction). The parser builds 46 throwaway instances per `nDMaterial` call. | minor / doc | none observed | reading; `test_H14_*` (structural) |

## §2 Adjudications where the lanes disagreed

- **H1 severity.** Red (both lanes) called it a results blocker; blue measured a
  two-cube heterogeneous model converging to each cube's stand-alone stress at 1e-6 with
  H1 fully active. **Blue is right on results: a converged step is exact because the
  residual is built from each element's own stress.** Red is right that it is a blocker
  for two other reasons: the 62 % / 5.3× iteration cost compounds with M3 on real meshes,
  and no threaded assembly (ADR-75b) can ever include this material. Filed as M1, not B.
- **DP NaN root cause.** R1-C blamed the `CHECK_APEX_REGION` stub and the dead apex call
  site; red-numerics showed the BE NaN guard fires and the host drops the −1, and traced
  the NaN to the `*= 0.0` pseudo-initialisation. **Both are true and neither alone is the
  fix**: the guard-drop is B2, the uninitialised storage is B4's one-line fix, the dead
  apex site is the reason a proper apex return cannot be reached even for HB, whose apex
  methods are implemented. Platform-dependent: Ubuntu CI (fresh heap) commits a clean
  finite history on the same path; only the dirty Windows pytest heap reproduces the
  NaN — which is the signature of UB, not of a deterministic defect.
- **`depsilon_elpl` blast radius.** R1-A said zero (no registered PF/YF reads `depsilon`);
  red-numerics found `ArmstrongFrederickPolicy` reads it and 22/46 specializations carry
  AF. **Red stands**; the static is stale or zero for every AF material under `Continuum`
  and `Secant`, and in the ME/RK45 drift correction. Folded into M3.
- **H6 on VonMises.** The plan's proposed order study was degenerate (R1-A); the tangent
  finding, not the accuracy claim, is the operational one. **D2 recommendation: do not
  rewrite in this ADR; document the map honestly and record that a closest-point BE is
  the only route to a consistent tangent.**
- **R4's Drucker–Prager and Hoek–Brown rows.** The harness's DP yield function was first
  written tension-positive, then negated after an empirical probe; `meanStress()` is
  `trace/3` (tension-positive, `OTHER/eigenAPI/typedefs.h:258`) so the source is
  self-consistent in OpenSees' convention and only the comment in `DruckerPrager_YF.h`
  is wrong. The agent's negation contradicts that reading, DP's |f| floor of ~1e-3
  relative is attributed to B3, and HB's 2.1e3 reflects M4. **The VM, MC and MCTC rows of
  `_adr94_matrix.md` are usable; the DP and HB rows are not warrant-grade and the
  "recommended configuration" list must not be read as accuracy guidance for those two.**
  The row that is usable: MCTC's best cell is `Backward_Euler` + `Continuum`, agreeing
  with ADR-84 §9.4.
- **R3a's VM "35 % derivative bug".** Withdrawn as a bug and re-filed as B5: VM is
  internally consistent under the tensor convention; the defect is the split across
  families and across consumption sites. DP's normal-slot error is independent and real.

## §3 Fail-loud contract state

| Integrator | Failure sites | Code | `TenNodeTetrahedron` | `LadrunoBrick` | `stdBrick` |
|---|---|---|---|---|---|
| Backward_Euler | strict exhaustion, strict SR-fallback | sentinel | propagated | propagated | swallowed |
| Backward_Euler | singular tangent, NaN, `dLambda<0` fallback | −1 / −1 / **0** | propagated / propagated / **none** | swallowed / swallowed / none | swallowed |
| FE, FE_sub, BE_LS, ME, RK45, RK45_old | all (13 sites) | −1 | propagated | **swallowed** | swallowed |
| all | f-decreasing elastic exit (8 sites) | 0 | none | none | none |

"Fail loud today" = `Backward_Euler` + `strict_convergence 1` + tet or `LadrunoBrick`
host, on 2 of ≥5 BE exit paths, in a unit system where the absolute tolerance binds.
Fix direction (R5): widen the sentinel to every site, do not loosen the hosts to `< 0`.

## §4 Conventions and bookkeeping

- Stress sign: tension-positive throughout (`meanStress = trace/3`); the DP comment
  saying otherwise is wrong. Shear slots: split (B5).
- Tolerances: absolute, stress-unit dependent (M5). Any ADR-94 follow-up must specify
  units next to every tolerance it quotes.
- Diagnostics: `cout`, not captured by pytest `capfd` on this build (child-process
  capture required; helper `_run_child` in `test_adr94_hlist_mechanical.py`).
- Ledger row 337 marker debt: markers added under this ADR (comment-only, same PR).
- No banner change (no feature ships). No class tag consumed.

## §5 Warrant status

- **Warrant-grade (runtime-pinned, order-independent, 3× rerun clean):** B1, B2, B3,
  B4 (NaN + guard), M1, M2 (4/6 sites), M3 (tangent), M4, M5, M6 (`reset`), M7, M8 (ME),
  m1.
- **Analytic-grade (exact by construction, runtime magnitude unmeasured):** B5.
- **Reading-only (not warrant-grade; red-tests §5):** H2, M2's BE_LS and RK45 sites,
  M6's `revertToLastCommit` half (tet self-heal hides it), M8's RK45 half, M9, m2, m3.
- **Sentinel gap (red-tests §1):** B1's tests gate only a control-flow fix; a warn-only
  fix would leave them green. Add a stderr assertion when fixing.
- **Coverage after this review:** runtime evidence exists for VM, DP, MC, MCTC, HB
  specializations on BE (and the R4 matrix on all six integrators for those five);
  StiffSoilCap has none, StiffSoilShear only NaNs; hardening-law combinatorics (Linear
  tensor, AF) untouched; `Numerical_Algorithmic_*` measured only on VM.

## §6 Scope and owner decisions

- **D1** (review only, fixes as separate WPs): held. No `SRC/` code changed on #804;
  comment-only markers excepted.
- **D2** (H6): recommend document + honest naming; no closest-point rewrite under this
  ADR. If a consistent tangent is ever wanted (M1 + M3 together are the 5.3× cost), open
  ADR 95.
- **D3** (SP/MP): out. M10 recorded only.
- **D4** (HB): port jaabell's `60d9b9b23`; ours is the older tree. Needs the exchange
  with José before any fork commit touches `HoekBrown_YF.h`.

## §7 Fix list

**(a) Fork-local, opt-in or default-inert, cheap — one WP each, sentinel tests flip:**

**Status:** the first four rows below SHIPPED together as **`wp/94a-fail-loud`** (one
coherent change: "the material can always fail loud, and misconfiguration cannot be
silent"), rather than as four separate WPs — they touch the same two files and the
sentinel widening is what makes the parser refusal and the strict-mode gates observable.
`wp/94e-revert` and the docs row are NOT done.

| WP | Fix | Flips |
|---|---|---|
| `wp/94a-parser-loud` **[SHIPPED in `wp/94a-fail-loud`]** | `opserr` + construction failure on unknown integration-option / model-parameter tokens; required-parameter assertion (ADR-84 P2(e)). Add a stderr assertion to the H13 tests. | `test_H13_*` |
| `wp/94b-sentinel-everywhere` **[SHIPPED in `wp/94a-fail-loud`]** | Return `LADRUNO_MATERIAL_REFUSED` at all 15 failure sites (incl. the `dLambda<0` fallback under `strict_convergence`, and the f-decreasing exits under the flag on every integrator). Keep hosts sentinel-only. | `test_H5_*`, `test_H7_*`, `test_R2_*stdbrick*` (expected still-swallowed on `stdBrick`), `test_adr94_contract` |
| `wp/94c-eigen-init` **[SHIPPED in `wp/94a-fail-loud`]** — six sites, not four (`DuncanChang_EL`'s `EE_MATRIX` and DP's degenerate `dev_part` branch too), plus 27 `*= 0` on ASDP's own class-statics | `setZero()` at the four `*= 0.0` / uninitialised `VoigtVector` sites (DP YF/PF, NullHardeningTensor, AF). Grep-gate the idiom. | `test_H10_dp_apex_*` |
| `wp/94d-refuse-bels` **[SHIPPED in `wp/94a-fail-loud`]** | Parser refusal of `Backward_Euler_LineSearch` and `Runge_Kutta_45_Error_Control_old`; drop dead enum values from the setter. | `test_H8_*` |
| `wp/94e-revert` **[NOT DONE]** | Implement `revertToLastCommit`/`revertToStart` (Trial ← Commit, IVs revert, `first_step` reset). | `test_H4_*`, contract revert tests |
| docs | `tangent_type Continuum` as the general recommendation; `strict_convergence 1` + BE required for softening; tolerance-in-units note. | — |

**(b) jaabell-bound framework work (fresh-branch ports, no AI traces, bit-identical gate;
coordinate before any fork default changes):**

| Item | Content |
|---|---|
| H1/F1/F2 | Per-instance `Stiffness`, `dsigma`, `depsilon_elpl`, `intersection_*`; `getInitialTangent()` without side effect; YF/PF `vv_out` by value or per-instance. One deliberate pass with a full-battery bit-identical check. |
| B5 | Pick one shear-slot convention for `df_dsigma_ij`/`pf` and one contraction at every site; fix DP's normal-slot gradient. Needs a numpy radial-return oracle on simple shear first (the runtime magnitude this review did not get). |
| M5 | Relative tolerance (`|Φ|` scaled by the YF's own strength scale). |
| B4 apex | Revive the apex call site + implement DP/RoundedMC apex methods (RoundedMC will not compile today). |
| M4 | HB composite (D4). **shipped wp/94d.** |
| M3/M8/M9 | Document the map; dead drift checks; per-stage `E`. Batch with H1. |

**(c) Documentation-only:** H2, H14/F6, m2's parameter-id map, M10.

## §8 Test results

`tests/`, `PYTHONPATH=../dist/bin python3.12 -m pytest` on the nine `test_adr94_*.py`
files plus the four pre-existing ASDP files: **69 passed, 4 skipped, 0 failed, 33.8 s**
(skips: 2 × `h5py` not installed in this shell, 2 × DP/HB/MCTC component smoke cells).
Baseline (R0, four pre-existing files alone): 22 passed. Red-tests lane: 27/27 in three
orderings and alone; three full reruns green. All sentinel tests assert the defect as
observed on `52314165a`; whoever fixes an item rewrites its test to assert the corrected
behaviour rather than deleting it.

## §9 Ledger obligations discharged in this PR

`LEDGER_quirks.md`: seven new entries (class-static tangent; sentinel-only host compare
drops bare −1; `capfd` vs `.pyd`; `WinError 6` under pytest; `printA` sparse-only;
tet `eleResponse` self-heal; absolute tolerance is unit-dependent).
`LEDGER_implementations.md`: one ADR-94 row (review, no class tag).
`LEDGER_vanilla_files.md`: row 337 updated (markers added). No new vanilla edits.
`94_asdplastic_review_plan.md`: implementation log + status flip. README line updated.
