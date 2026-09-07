# ADR-94 R2 — RED NUMERICS: the integrators/tangents are unsafe as shipped

Build `52314165a`. Reproducers `tests/test_adr94_redblue_numerics.py` (4 tests, 0.31 s;
26 passed with the three R1 ASDP files). Lines are `ASDPlasticMaterial3D.h` unless named.
See also `_adr94_hlist_R1A/R1B/R1C.md`, `_adr94_hb_drift.md`, plan §4.

**Thesis.** R1 showed checks are *missing*. The stronger fact: the checks that exist
**cannot reach the analysis** — on `stdBrick` no integrator can fail loud today,
including `Backward_Euler` with `strict_convergence 1`.

## Q1 — is the `strict_convergence` contract incomplete?

R1's eight accept-without-check sites: `compute_local_stress` 508 (`yf_start > yf_end`,
**no tol at all**), FE 1423, FE_sub 1599, BE_LS 2406, RK45_old 2667, ME 3088, RK45 3435,
and R1-A's eighth, BE's `dLambda + deltaLambda < 0 → return 0` (2298–2303). **Four more:**

- **9/10.** ME 3234 / RK45 3755 accept a substep unconditionally when
  `effective_dT <= dT_min`; 3253/3775 clamp `dT >= dT_min`, so a hard step *always*
  reaches it — the "error control" degenerates to fixed-step, accept-all.
- **11.** `return_to_yield_surface 1` (3270–3290, 3795, 2915, 1506, 1675) never verifies
  its result and silently no-ops when `abs(denominator) <= eps`.
- **12.** The two loud guards return a **bare `-1`** — `SINGULAR LOCAL TANGENT` 2280 and
  the NaN guard 2318–2323 — which every hex host drops.

**Can a user make ANY integrator fail loud? On `stdBrick`, no.** `Brick::update()`
assigns `success = mat->setTrialStrain(...)` then `return 0;` unconditionally, discarding
a bare `-1` **and** `LADRUNO_MATERIAL_REFUSED`. Measured, one deck/material/tolerance
(`test_R2_strict_convergence_is_a_noop_on_stdbrick`): `LadrunoBrick` **0/20**
(`analyze() == -3`), `stdBrick` **20/20, all `0`**. `LadrunoBrick` compares
`== LADRUNO_MATERIAL_REFUSED` in every formulation (std/SSP/URI/EAS/hypo/finite), so it
drops bare `-1`s too; only `TenNodeTetrahedron` (`success += …; return success`)
propagates both. Reachable set: **BE + `strict_convergence 1` + LadrunoBrick/tet**, on
2 of ≥5 BE exit paths.

## Q2 — `f_absolute_tol` is absolute: meaningful for HB/StiffSoil?

**No — the unit system decides pass/fail.** Default `1e-6`
(`OPS_AllASDPlasticMaterial3Ds.cpp:294`) on a `|Phi|` whose natural scale is `sigma_y`
(VM), `c cos(phi)` (MC), `sigma_ci*s^a` = 5350 kPa (HB at 50 MPa) — four orders across the
catalogue before any unit choice. Measured
(`test_R2_f_absolute_tol_makes_strict_convergence_unit_dependent`): the same MC problem in
**kPa** completes **20/20** at the default with `strict_convergence 1`; in **Pa** (×1000,
identical strains, identical physics) it is **refused on step 1, 0/20**. Tightening to
`1e-10` refuses *both* — the default sits in a ~4-decade window between "never binds" and
"never runs", and that window moves with the units.

For HB, `df_dsigma_ij` is a 6-point **central difference** with
`ds = max(HB_ds, HB_ds*||sigma||)` through a YF **discontinuous at `arg = 0`** (H10a):
inside a ball of radius `ds` the normal is `O(jump/2ds)` — `sigma_ci*s` = 587 kPa with
`ds` = 1e-4 gives `~3e6` of spurious gradient, worse the *finer* `HB_ds` is set. An
absolute 1e-6 is unreachable there by construction; the Newton then hits the H7 fallback
— R1-C's stall, mechanism now named.

## Q3 — is the BE internal-variable update consistent?

**No, except degenerately.** BE accumulates `trial_value += deltaLambda*h(sigma^k, iv^k)`
*inside* the Newton loop (2311–2316): `Delta_iv = sum_k dl_k h_k`, a quadrature over the
iterates, where true BE needs `Delta_lambda * h(sigma_{n+1}, iv_{n+1})`. Equal only if `h`
is constant over the iterates. `LinearHardeningForScalarPolicy` `h = H sqrt(2 m.m/3)` is
constant only for a *normalised* `m` (VonMises_PF) — DP/MC/StiffSoil PFs vary `||m||`.
`LinearHardeningForTensorPolicy` `h = H dev(m)` is inconsistent whenever the flow direction
rotates, i.e. every non-VM YF. `ArmstrongFrederickPolicy`'s `h` depends on `alpha` = the
IV's own pre-update trial value, so **AF's recovery term is integrated EXPLICITLY inside an
implicit return**, with a hard `if (alpha_norm >= alpha_limit) derivative = 0` clamp
(`AllASDHardeningFunctions.h:131-135`) that also makes `dPhi/dlambda` discontinuous.

So for everything but VM+linear the converged state is a function of the Newton *path*, not
of `(sigma_n, iv_n, Delta_eps)` — step-halving need not reproduce it. **22 of 46 registered
specializations carry AF.** Two more defects there, both missed: AF prints **6–7
unconditional `cout` per call** (per IV, per iteration, per GP, per step — H12 counted only
the commit line); and `eq_norm = sqrt(2/3 * v.squaredNorm())` uses Eigen's plain sum of
squares instead of `tensor_dot_stress_like`, under-weighting shear by 2.

## Q4 — `One_Step_Return`: consistent? does it update IVs?

**IVs yes; consistent no, three counts.** (1) **Wrong `Eelastic`**: in ME (3270–3290) it
uses the `Eelastic` left by the last substep's *corrector stage*, `et(predictor_sigma)` at
3178 — an intermediate, possibly rejected state, neither commit nor trial; same in RK45
(3795); in BE_LS `et(CommitStress)` is hoisted at 2379 outside the split loop. For
`StiffSoil_EL`/`DuncanChang_EL` the drift correction uses a discarded state's elastic
operator — H15 extended into the correction, which R1 did not examine.
(2) **`depsilon_elpl` is the class static**, passed to `pf(...)` *and*
`iv.hardening_function(...)`. R1-A called the blast radius zero because "no registered PF or
YF reads `depsilon`" — that grep excluded the hardening functions, and
**`ArmstrongFrederickPolicy` reads it**; BE never assigns `depsilon_elpl` (only 516/531,
1430/1441, 1606–1623, 2673/2684, 3095/3106, 3442/3453 do), so it is zero or another
integrator's leftover. This is a **stress-changing** path, not only a tangent.
(3) **One linearised step, never verified** — no loop, no post-check against `tol_yf`,
silent no-op if `abs(denominator) <= eps`. With H9's drift checks dead, the explicit
integrators have **no** verified admissibility mechanism at all. BE never calls it.

## Q5 — who else shares the dead apex call site (2093–2161)?

`yf_has_apex = true_type`: `DruckerPrager_YF`:137, `HoekBrown_YF`:220, `MohrCoulomb_YF`:219,
`RoundedMohrCoulomb_YF`:178, `TensionCutoff_YF`:149. **VM cannot** (no apex). **DP does.**
**MC does not** — measured (`test_R2_mc_hydrostatic_tension_through_apex_stays_finite`): 2×
the MC apex volumetric strain (`p_apex = c cot(phi)` = 173.2 kPa) commits finite bounded
stress, because MC's default (`MC_ds = 0`) analytical derivative guards `J2 < 1e-15`
explicitly; DP has no such guard. **HB** reaches its apex through the *discontinuous* branch
— it stalls rather than NaNs.

**R1-C's root cause is wrong.** R1-C: "nothing in the existing gates trips". Measured
(`test_R2_dp_apex_nan_guard_fires_and_is_swallowed`, child-process capture): the BE NaN
guard **does** fire — `"NaN!"` printed, `-1` returned — and `analyze()` still returns `0`
every step with NaN committed. **`LadrunoBrick` drops the `-1`.** H10b is a return-code
contract defect; fixing the apex would not fix the failure class.

**Where the NaN comes from (new).** `DruckerPrager_YF::df_dsigma_ij` builds
`VoigtVector pressure_part;` then `pressure_part *= 0.0; // Initialize to zero`
(`DruckerPrager_YF.h:64-65`; identically `DruckerPrager_PF.h:68-69`). `VoigtVector()`
forwards to Eigen's fixed-size default ctor, which does **not** zero, and
`EIGEN_INITIALIZE_MATRICES_BY_ZERO` is defined nowhere in this tree (grep) — so
`garbage*0.0` is NaN for a non-finite stale slot. Those are the **only two** `*= 0.0`
pseudo-initialisations in the whole tree, and they sit in the **only** YF that NaNs;
ADR-84 recorded this exact trap once already. Same UB twice more:
`NullHardeningTensorPolicy` returns an uninitialised `VoigtVector zero;`
(`AllASDHardeningFunctions.h:71-74`) added straight onto a back stress, and AF's
`VoigtVector derivative;` (:128) on the saturation branch. **Bonus:**
`RoundedMohrCoulomb_YF` declares the apex trait but defines **neither**
`check_apex_region` **nor** `apex_stress` — reviving the dead call site would not compile
for its registered specializations.

## Q6 — severity, Cerro-Lindo-scale MC/MCTC on `Backward_Euler` + `Secant`

| # | Finding | RESULTS or COST |
|---|---|---|
| 1 | **H1 class-static `Stiffness`** (R1-A): model assembled with the last GP's tangent | **RESULTS** (equilibrium path) + cost. Blocker |
| 2 | **`strict_convergence` cannot reach the analysis** (Q1) | **RESULTS**. Blocker; cheapest fix here |
| 3 | **H7** commits the elastic predictor with `f` growing, returns 0 | **RESULTS**. Blocker |
| 4 | **Tolerance is absolute** (Q2): kPa vs Pa decides refusal | **RESULTS** + cost. Major |
| 5 | **BE IV update is a Newton-path quadrature** (Q3) | **RESULTS**; MC is perfectly plastic so *this* config is spared, MCTC/AF/StiffSoil are not. Major, latent |
| 6 | **No tangent option is consistent; `Secant` default = 5.3× `Continuum`'s iterations** (H6) | **COST** only — stresses exact |
| 7 | **`depsilon_elpl` zero/stale in Continuum+Secant and in the ME/RK45 correction; AF reads it** (Q4.2) | cost today, **RESULTS** for the 22 AF specializations. Major |
| 8 | **DP `pressure_part` / `NullHardeningTensor` uninitialised-VoigtVector UB** (Q5) | **RESULTS** (NaN). Off the MC path; one line each. Major |
| 9 | **ME/RK45 accept-at-`dT_min`; `One_Step_Return` unverified with the wrong `Eelastic`** | **RESULTS**. Off-path for BE; blocks the explicit integrators. Major for StiffSoil |
| 10 | **BE_LineSearch ignores both options and truncates the step** (H8) | **RESULTS** if selected. Refuse the option |
| 11 | AF `cout` flood + `eq_norm` shear under-weighting (Q3) | cost / **RESULTS** under shear. Off-path |

## What R1 missed

1. **The gate is unreachable, not merely incomplete** — R1 never checked whether a present
   check survives the host. `stdBrick` 20/20 vs `LadrunoBrick` 0/20, same deck.
2. **R1-C's DP-NaN root cause is wrong** — the guard fires and is swallowed; the NaN's
   origin is a `*= 0.0` pseudo-init on an uninitialised `VoigtVector`, DP-only.
3. **R1-A's "`depsilon_elpl` blast radius is zero" is refuted** — the grep excluded
   `AllASDHardeningFunctions.h`; AF reads `depsilon`, and 22/46 materials carry AF.
4. **Tolerance semantics were never tested** — units alone flip 20/20 to 0/20.
5. **The BE IV update is not a BE IV update** for non-constant `h`.
6. **Four more accept-without-check sites** — ME/RK45's `dT_min` unconditional accept
   (with the clamp that guarantees you reach it), `One_Step_Return`, and the bare-`-1` guards.
7. **`RoundedMohrCoulomb_YF` declares an apex it does not implement** — a compile failure
   waiting for whoever revives the dead call site.
