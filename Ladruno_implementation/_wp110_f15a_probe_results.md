# WP-110 / F15a -- GetElastoPlasticTangent measurement probe: results

**Status: CONFIRMED.** `ManzariDafalias::GetElastoPlasticTangent`
(SRC/material/nD/UWmaterials/ManzariDafalias.cpp:5110) computes a wrong
continuum elastoplastic tangent at a plastic (loading) state. Mechanism:
the flow-direction vector `R` is converted contravariant -> covariant
**twice**:

```cpp
R = ToCovariant(temp0);                     // :5127 -- single, correct
temp1 = DoubleDot4_2(aC, ToCovariant(R));   // :5129 -- R is ALREADY covariant;
                                             //          this re-wraps it.
```

`ToCovariant` (:5651) only doubles the **shear** Voigt entries (indices
3,4,5 = xy,yz,xz); the normal entries (0,1,2 = xx,yy,zz) are untouched. Since
`aCep = aC - (Macauley/temp3) * Dyadic2_2(temp1, temp2)` scales its
correction term row-wise by `temp1`, the bug can only ever corrupt aCep's
**shear rows** (3,4,5 -- dSigma_xy, dSigma_yz, dSigma_xz w.r.t. any strain
component). The normal rows (0,1,2) are mathematically untouched. This is a
falsifiable, mechanism-specific prediction and it is exactly what was
measured -- see "Row-identity check" below.

Test: `tests/test_manzari_ep_tangent_probe.py` (passes; it records the
comparison, it does not gate on the bug being present/absent). Raw console
output: `Ladruno_implementation/_wp110_f15a_probe_raw_output.txt`.

## Build

```
ladrunoBuild = 48c0e99bc8e28bbb4fdf965015f285f01e16e90d
```

Matches `ladruno` tip `48c0e99bc` (>= the required `48c0e99`). This is the
**installed** Ladruno binary (`C:\Program Files\Ladruno\OpenSees\bin`), per
the WP-110 rule of no build in this worktree.

## Method (brief; full rationale in the test's module docstring)

- Single 8-node `SSPbrick` (1-point stabilized hex -> one shared material
  state) with `ManzariDafalias` in its native 3D form (no plane-strain
  reduction, so `eleResponse('stress')` exposes all 6 components including
  sigma_zz, which the 2D `PlaneStrain` wrapper hides).
- All 21 non-pinned DOFs driven by SP constraints under per-DOF `Path`
  TimeSeries (`constraints('Transformation')` -- `Plain` silently drops any
  non-homogeneous SP value, a real trap hit while building this harness) so
  the element sees an exact, fully prescribed affine (uniform) strain
  history with **zero free DOFs** -- no equilibrium solve can contaminate a
  probe.
- **FD leg:** central-difference the wrapper's 6-component stress against
  each of the 6 independent strain components, rebuilding the model from
  scratch per probe (replay the identical path, bump one extra checkpoint).
  This converges, as h -> 0, to the true *continuum* (rate-form) tangent
  regardless of which IntScheme produced the probe steps.
- **Analytic leg:** `GetStateDependent` / `GetElasticModuli` /
  `GetElastoPlasticTangent` transcribed verbatim into numpy (as-written and
  with the second `ToCovariant` removed), fed the exact internal state
  (mSigma, mAlpha, mFabric, mAlpha_in, void ratio, dGamma) read back from the
  same committed FE state via `eleResponse`. The 3D wrapper flips sign
  uniformly (`mSigma_M = -mSigma`); since both stress and strain flip
  together, the tangent matrix is invariant, so no sign correction is needed
  when comparing the two legs.
- IntScheme = 2 (`BackwardEuler_CPPM`), TanType = 1 (continuum aCep, not the
  separately-computed "consistent" tangent) throughout the plastic cases.

**Gotcha recorded (also in the test's docstring):** the full `ManzariDafalias`
constructor sets `mElastFlag = 0` ("stage 0") unconditionally, which forces
`elastic_integrator` regardless of IntScheme/TanType -- `GetElastoPlasticTangent`
is never reached and TanType 0/1/2 are silently identical until
`ops.updateMaterialStage('-material', tag, '-stage', 1)` is called. Cost real
time to find; the test asserts `dGamma > 0` at the plastic states specifically
to catch a regression back into this trap.

## Elastic sanity gate

TanType=0, small isotropic strain, no plasticity possible:

| | analytic Ce | FD |
|---|---|---|
| normal diag | 7198 | 7198 |
| normal off-diag | 378.8 | 378.8 |
| shear diag | 3410 | 3410 |

Max relative error: **2.2e-9**. Confirms the FE harness and the numpy
transcription of `GetStiffness`/`GetElasticModuli` are both correct before
touching the suspect function.

## Plastic-state results

Two states reached via a 2-stage strain ramp (isotropic consolidation, then
an added deviatoric/shear increment), IntScheme=2, TanType=1:

| | LOW p' | HIGH p' |
|---|---|---|
| p' (kPa) | 205.9 | 13218 |
| q (kPa) | 178.4 | 10976 |
| eta = q/p' | 0.866 | 0.830 |
| dGamma | 6.48e-5 | 9.53e-4 |

(dGamma > 0 confirms an active, loading plastic state at both.)

### Row-identity check (the falsifiable prediction)

At both states: **rows 0-2 (normal) of aCep are bit-identical between the
as-written (buggy) and fixed formulas; rows 3-5 (shear) differ.** This is
exactly the row-selective corruption pattern the double-`ToCovariant`
mechanism predicts, and rules out an unrelated/general-purpose bug.

Off-diagonal shear-row entries (row 3 or 4 or 5, column in 0-2) come out at
**exactly 2x** between buggy and fixed (e.g. LOW p', row 3: buggy
[27560, -20770, 28970] vs fixed [13780, -10380, 14480] -- ratio 2.000, 2.001,
2.001) — the aC contribution to those entries is zero (`GetStiffness` has no
normal-shear coupling), so they come ONLY from the rank-1 correction term,
which scales linearly with the doubled `temp1` shear component.

### Per-block relative error, FD vs analytic (h=1e-6, not yet fully h-converged)

| block | LOW p': FD vs buggy | FD vs fixed | HIGH p': FD vs buggy | FD vs fixed |
|---|---|---|---|---|
| normal-normal (rows/cols 0-2) | 7.6%-75% | *(identical to buggy — unaffected rows)* | 1.0%-99.8% | *(identical)* |
| shear-shear diag, e.g. (3,3) | 19.1% | 4.3% | 30.3% | 7.1% |

Rows 0-2 are identical between the "FD vs buggy" and "FD vs fixed" tables by
construction (those rows are bug-free in both variants); their nonzero
residual against FD reflects the known theoretical gap between the
*continuum* (rate-form) tangent `GetElastoPlasticTangent` computes and the
*algorithmic* tangent that a finite difference of the actual (discretely
integrated) response measures -- expected, not further bug evidence.

### FD step-size convergence, entry (3,3) = dSigma_xy/dEps_xy

LOW p':

| h | FD | buggy | fixed |
|---|---|---|---|
| 1e-4 | 37293 | 42635 | 48676 |
| 1e-5 | 46331 | 42635 | 48676 |
| 1e-6 | 50787 | 42635 | 48676 |
| 1e-7 | 51269 | 42635 | 48676 |
| 1e-8 | 51318 | 42635 | 48676 |

HIGH p':

| h | FD | buggy | fixed |
|---|---|---|---|
| 1e-4 | 338400 | 305800 | 372110 |
| 1e-5 | 391690 | 305800 | 372110 |
| 1e-6 | 398590 | 305800 | 372110 |
| 1e-7 | 399240 | 305800 | 372110 |
| 1e-8 | 398920 | 305800 | 372110 |

FD converges cleanly and monotonically to a stable plateau at both states
(the analytic values, being h-independent by construction, are flat lines).
At the converged limit:

- LOW p': **fixed is 5.1% off FD; buggy is 16.9% off** (buggy error ~3.3x fixed's).
- HIGH p': **fixed is 6.8% off FD; buggy is 23.4% off** (buggy error ~3.4x fixed's).

Fixed is consistently, substantially closer to the converged FD value than
buggy, at both confinement levels, for exactly the row family the mechanism
predicts.

## Verdict

**CONFIRMED.** Mechanism = double `ToCovariant` on the flow-direction term
`R` at ManzariDafalias.cpp:5129 (the loading-direction conversion at ~:5132
is single and correct, as suspected). Effect: aCep's shear rows (dSigma_xy,
dSigma_yz, dSigma_xz w.r.t. any strain component) get a rank-1 plastic
correction inflated by exactly 2x relative to the correct value; normal rows
(dSigma_xx, dSigma_yy, dSigma_zz) are completely unaffected. Measured
magnitude at the two states probed here is larger than the ~9% prior figure
(here: buggy is 17%-23% off the converged FD limit vs fixed's 5%-7%), and the
error appears state-dependent (grows somewhat with p'/eta in this sample) --
consistent with a real, non-trivial tangent defect that a Newton solver using
TanType=1 on any IntScheme calling `GetElastoPlasticTangent` (BackwardEuler
CPPM, RungeKutta variants, MaxStrainInc/MaxEnergyInc) would feel as a
systematically wrong, direction-dependent (shear-only) Jacobian.
