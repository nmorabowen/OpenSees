# ADR-94 R3a — component derivative harness (results)

Phase R3a of the `ASDPlasticMaterial3D` review (`94_asdplastic_review_plan.md`
Sec 5 "R3", first three bullets). HEAD `52314165a`, worktree
`asdplastic-review-plan-62585c`.

## Route taken

**Standalone C++ route (preferred route) worked end-to-end** — no fallback
needed. The YieldFunction / PlasticFlow / Elasticity headers under
`SRC/material/nD/ASDPlasticMaterial3D/{YieldFunctions,PlasticFlowDirections,
ElasticityModels}/` are header-only templates over Eigen with no OpenSees
core dependency actually exercised by `operator()`/`df_dsigma_ij` (the only
place `Vector`/`Matrix`/`OPS_Stream` leak in is `EigenAPI.h`'s
`VoigtVector::fromStrain/toStrain`-style conversion helpers, which the
harness never calls, and `Vector::operator()` is `inline` in `Vector.h`
itself, so no `.lib`/`.cpp` link is required). Compiled and linked cleanly
with MinGW g++ 15.2.0, `-std=c++17`, against:
  - Eigen from the conan cache used by the CMake build:
    `C:/Users/nmb/.conan2/p/eigen5481853932f72/p/include/eigen3` (confirmed
    from `build/build/Release/generators/Eigen3-release-x86_64-data.cmake`).
  - Every directory under `SRC/` and `OTHER/` that contains a `.h` (same set
    CMake's `file(GLOB_RECURSE ...)` wires in — see root `CMakeLists.txt`
    lines ~403-420), captured once in `adr94_oracle/incdirs.txt` /
    `incflags.txt`.
  - `using namespace std;` added locally in the harness (upstream relies on
    `StiffSoil_HardeningFunctions.h`'s stray file-scope `using namespace
    std;` leaking into the real TU's include order for unqualified
    `cout`/`endl`/`tuple` in several headers — see Trap below).

No `SRC/` edits, no CMake changes, no commits.

**Dead-code note (not a defect, but worth recording):** the two existing
"standalone harness" files named in the task brief,
`YieldFunctions/test_00_VonMises_YF.cpp` and
`PlasticFlowDirections/test_pf.cpp`, do **not** compile against the current
headers — they reference `VonMisesRadiusIV<...>` / `BackStressIV<...>`,
which no longer exist (today's names are `YieldStress<HardeningType>` /
`BackStress<HardeningType>`, see `AllASDInternalVariableTypes.h`). They are
stale leftovers from an earlier refactor, not wired into any build, and were
not usable as-is. This is the same class of finding already flagged for
`test_HoekBrown.cpp` in the review plan Sec 5 R3 (4th bullet).

## Harness

`Ladruno_implementation/adr94_oracle/fd_components.cpp` (build command in its
header comment). For each of the 7 registered YF families it instantiates
the **exact template combo used by the registry's first matching
registration** (copied from `ASD_material_definitions.cpp`), builds a
194-point stress cloud (160 random general 6-component states + 24
Lode-edge points at theta = ±30° for p ∈ {1,5,20,50} × r ∈ {2,10,30} + 6
hydrostatic-axis points p ∈ {0,1,5,20,50,100} + 4 J2→0 near-hydrostatic
points), and compares `df_dsigma_ij` against a central finite difference of
`operator()` (h = 1e-6), reporting the max relative error **split by
category**: "smooth" (random + Lode-edge, where the surface is genuinely
differentiable) vs. "singular" (the exact hydrostatic axis and J2→0, where a
`sqrt(J2)`-type surface has a real cone corner and no single gradient can
match a central difference approaching from an arbitrary direction — this is
expected non-smoothness, not a bug, and is reported separately so it does not
mask or get confused with genuine smooth-region errors). PF checks `m`
finite everywhere, and for VM/DP (etabar=eta) that `m` equals the YF's own
normal direction. EL checks `E(sigma)` symmetric + SPD (min eigenvalue) over
the cloud (StiffSoil_EL restricted to its admissible p > 0 subset).

## Results (smooth-region — the meaningful number)

Regenerated on `wp/94c-numerics` (commit `3622d6214`) with the same harness, same
194-point cloud, same seed. "vs Voigt FD" is the number that matters: after wp/94c
every yield function is differentiated in ONE convention (Voigt / engineering shear),
so the whole column is at finite-difference truncation noise. The "vs tensor FD"
column is the same gradient scored against the OTHER convention, and is now uniformly
~0.7 — i.e. the catalogue is no longer split, it is uniformly Voigt.

| Component | vs Voigt FD, pre-94c | vs Voigt FD, **wp/94c** | vs tensor FD, wp/94c | singular region (expected) | non-finite? |
|---|---:|---:|---:|---:|---|
| VonMises_YF | **3.53e-01** | **1.13e-08** | 7.06e-01 | 5.68e+00 | no |
| DruckerPrager_YF | **9.71e-01** | **1.15e-08** | 6.85e-01 | 2.97e+00 | no |
| MohrCoulomb_YF | 3.41e-06 | 3.41e-06 | 6.54e-01 | 2.58e-01 | no |
| HoekBrown_YF | 2.45e-06 | 2.45e-06 | 7.06e-01 | 1.57e-01 | no |
| StiffSoilCap_YF | 2.18e-06 | 2.18e-06 | 6.95e-01 | 9.83e+01 | no |
| StiffSoilShear_YF | 1.52e-05 | 1.52e-05 | 7.06e-01 | 9.90e-01 | no |
| MohrCoulombTensionCutoff_YF | 2.15e-07 | 2.15e-07 | 6.54e-01 | 8.91e-01 | no |

VonMises and DruckerPrager are the only two rows that moved, and they are the only two
wp/94c touched. VM's 3.53e-01 was the *convention* split (it returned the bare tensor
derivative); DP's 9.71e-01 was a genuine gradient error on top of it — `d sqrt(J2)/d v`
is `r/(2 sqrt(J2))` on the normal slots and `r/sqrt(J2)` on the shear slots, and the
code applied the shear answer to all six. The five analytical-or-numerical Voigt YFs are
bit-for-bit unchanged, which is the evidence that wp/94c moved VM/DP *onto* the majority
convention rather than moving the goalposts.

The two dropping to ~1e-8 rather than ~1e-6 is not a better derivative: VM and DP are
now CLOSED-FORM against a closed-form FD, while MC/HB/MCTC/StiffSoil differentiate
their own numerically differentiated `df_dsigma_ij`, so their floor is FD-on-FD noise.

PF / EL:

| Component | Result |
|---|---|
| VonMises_PF (associated, alpha shared with YF) | 0/194 non-finite; **max_assoc_err = 0.0** (exactly matches YF normal) |
| DruckerPrager_PF (etabar = eta, associated) | 0/194 non-finite; **max_assoc_err = 0.0** (exactly matches YF normal — pre-94c it inherited the SAME 2x-normal-slot defect as the YF, consistently; wp/94c fixed both, and the associativity identity still holds exactly) |
| MohrCoulomb_PF | 0/194 non-finite |
| HoekBrown_PF | 0/194 non-finite |
| StiffSoilCap_PF | 0/194 non-finite |
| StiffSoilShear_PF | **6/194 non-finite** (see below) |
| MohrCoulombTensionCutoff_PF | 0/194 non-finite |
| LinearIsotropic3D_EL | symmetric to machine precision; min eigenvalue 1.15e4 (SPD, 0/194) |
| StiffSoil_EL | symmetric to machine precision; min eigenvalue 1.11e4 (SPD, 0/113 admissible points) |

## Findings

1. **REVISED (see "Convention adjudication" below) -- NOT a gradient bug** — `VonMises_YF::df_dsigma_ij` under-weights shear
   components by a factor of 2.** `operator()` computes
   `f = sqrt(Q) - ...` with `Q = tensor_dot_stress_like(dev,dev) = sum(normal^2) + 2*sum(shear^2)`
   (the "stress-like" dot product doubles shear cross-terms to recover a true
   tensor double-contraction — see `EigenAPI.h`/`typedefs.h`). The true
   gradient is `dQ/d(shear_i) = 4*shear_i` vs. `dQ/d(normal_i) = 2*normal_i`,
   so `d(sqrt(Q))/d(shear_i) = 2*shear_i/sqrt(Q)` vs.
   `d(sqrt(Q))/d(normal_i) = normal_i/sqrt(Q)`. The code
   (`VonMises_YF.h` ~line 68) returns `dev/den` uniformly for all six
   components — correct for the normal components, **half the correct value
   for the three shear components**. Measured smooth-region relative error up
   to 35% at shear-heavy stress states; exactly zero at pure-normal states
   (consistent with the mechanism). `VonMises_PF` reuses the identical
   (buggy) formula, so the associated-flow equality check trivially passes
   (both sides are wrong in the same way) — associativity is preserved, but
   both the yield normal and the flow direction are wrong off the pure-normal
   axis.
2. **CONFIRMED defect, convention-independent (see adjudication below)** — `DruckerPrager_YF::df_dsigma_ij` over-weights normal
   components by a factor of 2 (mirror image of #1).** Here `sqrt_J2 =
   sqrt(0.5*Q)`, so the true gradient is `d(sqrt_J2)/d(normal_i) =
   normal_i/(2*sqrt_J2)` vs. `d(sqrt_J2)/d(shear_i) = shear_i/sqrt_J2`. The
   code again returns `dev/den` uniformly — this time **correct for shear,
   2x too large for normal components**. At a pure-normal Lode-edge point
   (zero shear) this is a clean, exact 2x error: measured relative error
   0.971 (→ 1.0 in the limit). `DruckerPrager_PF` (with `etabar = eta`,
   the associated-flow configuration) reuses the same formula, so again
   `m == df_dsigma_ij` by construction — the pair is internally consistent
   but both sides are wrong for a general (non-pure-normal) stress state.
   These two findings should be escalated to the H-list / red-blue lanes
   (R2/R6) as a single root cause: **both YFs apply a flat `dev/den`
   normalization to a `sqrt(alpha * tensor_dot_stress_like(dev,dev))`
   expression without accounting for how `alpha` (1 for VM, 0.5 for DP)
   changes the per-component chain-rule factor between the "doubled" shear
   and "undoubled" normal Voigt slots.** MC/HB/StiffSoil/MCTC do not have
   this defect (all smooth-region errors ~1e-6..1e-7, consistent with plain
   FD truncation noise at h=1e-6) — they compute their gradients through a
   different (non `dev/den`-only) path per component.
3. **`StiffSoilShear_PF` returns non-finite for 6/194 cloud points.** Root-caused
   in ADR-97 P5 (`Ladruno_implementation/reviews/adr97_p5_report.md`): the 6
   points are exactly the six **hydrostatic-axis** points in the cloud
   (`sigma = (p,p,p,0,0,0)` for `p` in `{0,1,5,20,50,100}`) — NOT "random
   general (non-diagonal) points" as this bullet originally said.  At an
   exactly hydrostatic trial stress `computeMobilizedDilatancy()` returns
   `psi_m = 0` exactly, which collapses the central-difference numerator of
   `PLASTIC_FLOW_DIRECTION`'s flow-direction gradient to the exact zero vector
   on all six Voigt axes, and the unguarded `vv_out /= norm` at `norm == 0.0`
   is the indeterminate `0/0`. Fixed by a zero-guard on the normalization
   (`StiffSoilShear_PF.h`); the cloud now reads 0/194 non-finite.
4. **Singular-axis numbers are expected, not defects.** Every `sqrt(J2)`-type
   YF shows a large "singular" column because a central difference straddling
   the exact cone tip of the yield surface cannot match any single
   sub-gradient — this is a property of the surface, not the code. It is
   reported separately specifically so it is not mistaken for evidence
   against MC/HB/StiffSoil/MCTC's otherwise-clean smooth-region numbers.
5. `StiffSoilCap_YF`'s huge singular-column number (98.3) is driven by
   evaluating at the literal `p=0` corner where its friction/cap terms are
   most degenerate; same expected-non-smoothness caveat as #4, not
   independently investigated further here.
6. Elasticity: both registered `EL` models (`LinearIsotropic3D_EL`,
   `StiffSoil_EL`) are symmetric to machine precision and SPD (positive
   minimum eigenvalue, ~1.1e4) over their respective admissible stress
   ranges — no defect found.

## CI pin (openseespy route, `tests/test_adr94_components.py`)

Per the task brief, a Zone-A pytest re-observes the SAME signal without
depending on the standalone `.cpp` (which needs a hand-wired g++ + conan
Eigen path, not run in CI): build a free-DOF `stdBrick` unit cube, push it to
a converged plastic state, and compare the assembled tangent under
`tangent_type Continuum` (uses `df_dsigma_ij` inside the consistent-tangent
formula) against `Numerical_Algorithmic_FirstOrder` (numerically
differentiates the material's own stress response, so it never touches
`df_dsigma_ij` and cannot inherit its bug).

- **VonMises** (hard pin): converges at `load_z=-20, 20 steps`; observed
  Continuum-vs-NumAlgFirstOrder relative tangent mismatch **1.4%-1.8%**
  (smaller than the 35% seen directly on `df_dsigma_ij` because this uniaxial
  loading path is mostly normal-deviatoric, diluting the shear-only defect).
  Pinned `> 0.005`.
- **MohrCoulomb** (hard pin): converges cleanly; mismatch is **exactly
  0.0**, consistent with the harness's ~1e-6 (FD-noise-level) smooth-region
  error. Pinned `< 0.005`.
- **DruckerPrager, HoekBrown, MohrCoulombTensionCutoff**: exercised (smoke —
  no crash / hang) but **not threshold-asserted**. The admissible load window
  between "still elastic" (mismatch trivially 0) and "Backward_Euler fails
  to converge" proved narrow for the parameter sets tried
  (`DP_xi_c=5, DP_eta=0.3`; `HB_sigci=30, HB_mb=2, HB_s=0.01, HB_a=0.5`;
  `MC_phi=30, MC_c=10, TC_min_stress=-5`) within R3a's time-box; DruckerPrager
  in particular converged at `load_z=-2.0` (mismatch 0.0, still elastic) and
  at `load_z=-2.5` (mismatch 0.0) but failed to converge at `-3.0` and above.
  Their authoritative FD numbers are the ones in the table above, produced
  directly against `df_dsigma_ij` by the standalone harness (which has no
  such convergence constraint since it never runs a Newton loop). Widening
  the DP/HB/MCTC load window to also pin them through openseespy is left as
  follow-up.

Test result: **3 passed, 2 skipped**, wall time < 1 s (`pytest
tests/test_adr94_components.py -v`).

## Files

- `Ladruno_implementation/adr94_oracle/fd_components.cpp` — the harness
  source (compiles clean, ran successfully; not wired into CI).
- `Ladruno_implementation/adr94_oracle/incdirs.txt` / `incflags.txt` — the
  captured `-I` list (mirrors CMake's own glob) for reuse by any future
  standalone ASDPlasticMaterial3D header harness.
- `tests/test_adr94_components.py` — the Zone-A CI pin described above.

## Convention adjudication (coordinator-requested followup)

`tensor_dot_stress_like`/`tensor_dot_strain_like` (`typedefs.h`) both double
shear cross-terms (`v12*v12*2 + ...`) to recover a true tensor double
contraction from raw (undoubled) Voigt storage; `tensor_dot_energy_like` is
a plain `.dot()`. A scalar function of a symmetric tensor has two legitimate
per-component shear derivatives: **tensor** (∂f/∂σ12, unambiguous, no
symmetric-partner doubling) and **Voigt** (∂f/∂v12 = 2·tensor, since v12
drives both σ12 and σ21). My central difference perturbs the single stored
`v12`, so it measures the **Voigt** derivative unconditionally — this part of
the original write-up was correct as a measurement, but "off by 2 = bug" was
not adjudicated against consumption. Diagonal (normal) slots have NO such
ambiguity (tensor = Voigt always) — any normal-component mismatch is a real,
convention-independent bug.

**Second FD mode** added to `fd_components.cpp` (`max_rel_err_smooth_tensor`
= analytic vs. Voigt-FD-with-shear-halved): VonMises_YF now matches the
**tensor** convention to 1.3e-8 (FD noise) — its "35% bug" is a legitimate
convention choice, not an error. DruckerPrager_YF matches **neither**
(0.971 vs both) — its normal-slot error is real and convention-independent,
confirmed. MohrCoulomb/HoekBrown/StiffSoilCap/StiffSoilShear/MCTC match
**Voigt** to ~1e-6 and tensor to ~0.65-0.71 — they use the Voigt convention.

**Per-site consumption table** (`ASDPlasticMaterial3D.h`, grepped every
`n`/`m` contraction):

| Site (function, lines) | Contraction on `n^T·E·m` | Convention it requires for `n` |
|---|---|---|
| `ComputeTangentStiffness` (438-482) L462,476 | `n.transpose()*Eelastic*m` (plain) | Voigt |
| `compute_local_stress` (483-585) L546 | `n.transpose()*Eelastic*m` (plain) | Voigt |
| `Forward_Euler` (1386-1561) L1454 | `n.transpose()*Eelastic*depsilon_elpl` (plain) | Voigt |
| `Forward_Euler_Subincrement` (1562-1721) L1636 | `n.transpose()*Eelastic*m` (plain) | Voigt |
| `Backward_Euler` (2034-2354) L2275 | `tensor_dot_stress_like(n, Eelastic*m)` (doubled) | **Tensor** (old plain form is commented out immediately above, L2274) |
| `Backward_Euler_LineSearch` (2359-2631) L2606 | `n.transpose()*Eelastic*m` (plain) | Voigt |
| all integrators: `TrialPlastic_Strain += dLambda*m`, `TrialStress -= dLambda*Eelastic*m` | plain accumulation / matrix-vector | `m` must be Voigt/engineering (`Eelastic`'s shear diagonal is `mu` not `2*mu` — confirmed `LinearIsotropic3D_EL.h:59-61` — so it expects an engineering-shear strain input; OpenSees' own strain Voigt convention is engineering, per `hex8_tangent.py`'s docstring) |

**Verdict, revised per family:**
- **VonMises**: `df_dsigma_ij`/`VonMises_PF` are self-consistently **Tensor**-convention — not a gradient bug. But this is REAL and consequential in two ways: (1) 5 of 6 sites above (everything except the current `Backward_Euler`) pair `n` against `Eelastic*m` with a plain contraction that only recovers the correct scalar for Voigt-convention `n` — for VM's Tensor `n` they under-count shear's contribution to the plastic modulus/consistent tangent (only `Backward_Euler` L2275 is paired correctly for VM); (2) `m` is consumed everywhere as a direct strain-Voigt accumulator/matrix operand expecting Voigt/engineering convention, but VM's `m` is Tensor-convention — under-counts the plastic shear strain increment and its stress correction by 2x in **every** integrator, VM only.
- **MohrCoulomb/HoekBrown/StiffSoilCap/StiffSoilShear/MohrCoulombTensionCutoff**: Voigt-convention `n`/`m` — correctly consumed by the 5 plain-dot sites and by the `m`-accumulation logic. **The current `Backward_Euler`'s L2275 doubled contraction is therefore backwards for this family** — it over-counts shear's contribution to the plastic modulus for MC/HB/StiffSoil/MCTC, the opposite failure mode from VM. This one line cannot be correct for both conventions at once.
- **DruckerPrager**: independent of the above — a genuine, convention-blind 2x error confined to the normal (diagonal) `df_dsigma_ij` components (matches neither Tensor nor Voigt FD). Stands as reported.

**Runtime probe**: attempted a simple-shear `stdBrick`+VM cube (lateral load,
zero normal load, `Backward_Euler`+`Continuum`) per the coordinator's ask;
the rig did not reach a clean converged plastic state within the time
available (Newton stalled partway into the load ramp on this loading path —
a separate, unexplained convergence issue, not investigated further here).
The analytic FD-vs-both-conventions comparison above is exact and does not
depend on this probe; it is the decisive evidence.
