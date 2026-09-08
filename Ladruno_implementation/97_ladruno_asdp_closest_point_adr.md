---
title: ADR 97 — ASDPlasticMaterial3D closest-point return map with a consistent tangent
project: Ladruno
status: draft
priority: high
owner: nmora
tags:
  - implementation
  - material
---

# ADR 97 — ASDPlasticMaterial3D: closest-point return map with a consistent tangent

## What

A new, opt-in constitutive integrator for `ASDPlasticMaterial3D` —
`integration_method Closest_Point` — that solves the fully implicit return map

    sigma_{n+1} = sigma_trial - dlambda * E : m(sigma_{n+1}, q_{n+1})
    q_{n+1}     = q_n + dlambda * h(sigma_{n+1}, q_{n+1}, m)
    f(sigma_{n+1}, q_{n+1}) = 0

as one coupled Newton system in `x = (sigma_{n+1}, q_{n+1}, dlambda)`, plus the
exact consistent (algorithmic) tangent of *that* map, exposed as the previously
dead `tangent_type Algorithmic`. Scope is family-by-family: P1 the smooth
families (VonMises, Drucker–Prager including the apex) with Linear/AF/Null
hardening; P2 the principal-stress-space multi-surface families (MohrCoulomb,
MohrCoulombTensionCutoff); P3 HoekBrown. Out of scope: StiffSoil (two-surface)
and RoundedMohrCoulomb, which have no runtime coverage at all (ADR-94 §5) and are
gated behind P5; any change to the *default* integrator or tangent (a separate PR
gated on P6); MP/database `sendSelf` (ADR-94 D3, still out).

## Why

ADR-94's verdict ([[reviews/adr94_verdict]] §1 M2/M3/M9, §6 D2) measured the
gap and explicitly deferred the rewrite to this ADR:

- **M3 — no tangent is the tangent of the map.** Against a central difference of
  the material's own committed response: `Continuum` 57 %, `Secant` (the
  DEFAULT) 80 %, `Elastic` 103 %, `Numerical_Algorithmic_*` 31 % — the last one
  because it differentiates `compute_local_stress()`, a *third* map that nothing
  commits. The default costs 5.3x `Continuum`'s global iterations.
- **M3 — the shipped `Backward_Euler` is an Ortiz–Simo cutting plane**, not a
  closest-point map: it re-evaluates `n`, `m`, `H` at the running trial state and
  accumulates `sigma -= deltaLambda * (E*m)` per iteration, so its fixed point is
  `sigma_tr - sum_k dl_k E m(sigma^k)`, not `sigma_tr - dl E m(sigma_{n+1})` —
  the two coincide only when the flow direction does not rotate. Its IV update is
  a Newton-path quadrature, exact only for constant `h`; **AF's recovery term is
  integrated explicitly inside the implicit return** (22 of 46 registered
  specializations carry AF).
- **M2 — explicit-integrator drift.** FE, FE_sub, ME, RK45 commit `f_MC` in the
  hundreds-to-thousands; their drift checks are empty `if` bodies (M8); ME/RK45
  accept unconditionally at `dT_min`.
- **M9 — stress-dependent elasticity is integrator-dependent.** BE evaluates
  `E(sigma)` once, at commit; ME/RK45 per stage; `One_Step_Return` uses a
  discarded stage's `E`. StiffSoil/DuncanChang are different materials per
  integrator.

M1 + M3 together were the 5.3x/62 % Newton cost on the ADR-94 two-cube model;
M1 is fixed (#815), M3 is not.

## Where

- **Modify** `SRC/material/nD/ASDPlasticMaterial3D/ASDPlasticMaterial3D.h`: new
  `Closest_Point(const VoigtVector&)` member; a dispatch case at the
  `setTrialStrainIncr` switch (~line 292); an `Algorithmic` branch in
  `ComputeTangentStiffness()` (~line 477); a `Closest_Point` gate on
  `Algorithmic`; the `Numerical_Algorithmic_*` re-point (P4).
- **Modify** `YieldFunctionBase.h`, `PlasticFlowBase.h`, `HardeningFunction.h`,
  `ElasticityBase.h`: new macro-declared interface members (below), each with an
  inert base default plus an opt-in trait, so an un-converted family compiles
  unchanged and is *refused at parse time*, never silently approximated.
- **Modify** the per-family headers `YieldFunctions/{VonMises,DruckerPrager,
  MohrCoulomb,MohrCoulombTensionCutoff,HoekBrown}_YF.h`, the matching
  `PlasticFlowDirections/*_PF.h`, and `AllASDHardeningFunctions.h`.
- **Modify** `OPS_AllASDPlasticMaterial3Ds.cpp`: `integration_method` gains
  `Closest_Point`; `tangent_type` gains `Algorithmic`; a new
  `experimental_integrator` (int, default 0) option; the help string and
  `ASDP_VALID_INTEGRATION_OPTIONS`.
- **No new file, no new class, no new class tag.** ASDP is one templated class;
  `Closest_Point` is a member, exactly as `Backward_Euler` is. No new CMake
  target (the ASDP TU already needs `/bigobj`).
- **Marker rule:** every hunk carries `// Ladruno (ADR-97) ...` so the vanilla
  footprint stays reconstructable by `grep -rn "Ladruno" SRC/`
  (`LEDGER_vanilla_files.md` row per file, same PR).
- **Oracles:** `Ladruno_implementation/adr97_oracle/` (numpy/g++, no OpenSees
  link — the `adr94_oracle/` pattern). **Tests:** `tests/test_adr97_<phase>.py`.

## Decisions

**D1 — `Closest_Point` is a NEW opt-in method; `Backward_Euler` stays
byte-identical.** Every result the fork has published on this material (Cerro
Lindo M-series, the ADR-84 battery, the finite-strain D3 oracle) was obtained on
`Backward_Euler`. ADR-84 §3's bit-identical strategy is what keeps those results
citable after five fix PRs, and it is cheap here: the new map is a new member
function and a new switch case, touching no line the old one reads. The default
flip (`Closest_Point` + `Algorithmic` as shipped defaults) is a separate PR,
opened only if P6 measures a net win at mesh scale.

**D2 — `tangent_type Algorithmic` is the exact consistent tangent of the
`Closest_Point` map, and is refused with any other integrator.** The enum value
already exists and has no dispatch case anywhere — it is dead, which is why
nobody has been silently getting it. A consistent tangent is defined *only*
relative to a specific committed map; offering `Algorithmic` on the cutting-plane
`Backward_Euler` would ship a fourth almost-right tangent, which is exactly the
class of defect M3 found. The parser refuses `tangent_type Algorithmic` unless
`integration_method Closest_Point`, naming both tokens. Symmetrically, P4
re-points `Numerical_Algorithmic_*` to central-difference the **actual committed
map** instead of `compute_local_stress()` — closing ADR-84 P4 and M3's
"differentiates a third map", and making `Numerical_Algorithmic_FirstOrder` a
cross-check on `Algorithmic` rather than a competitor.

**D3 — family-by-family return maps, three algorithm classes.** The smooth
families (VonMises, Drucker–Prager, including the apex) get the standard CPPM in
6D Voigt (Simo & Hughes ch. 3; de Souza Neto, Peric & Owen 2008 ch. 7–8), where
`f` is C1 away from the vertex and the Jacobian blocks are analytic. MohrCoulomb
and MohrCoulombTensionCutoff get a principal-stress-space multi-surface CPPM with
face / edge / apex (/ cutoff-corner) returns and the Koiter consistent tangent
(Clausen, Damkilde & Andersen 2006, 2007), reusing ADR-84's `special_return`
geometry and its SR tangent slot. HoekBrown gets the principal-space return of
Clausen & Damkilde (2008). The split is forced by the surfaces: MC's 6D `f` is
the exact Lode-angle form `A(theta)*sqrt(J2) + sin(phi)*p - c*cos(phi)`, whose
gradient carries a `1/cos(3*theta)` that the shipped code dodges with a
`|theta| < 29 deg` guard and a Drucker–Prager smoothing branch — a Newton system
built on that gradient cannot converge quadratically at the corner, which is
precisely where MC models spend their time; in principal space those corners are
exact linear features. StiffSoil (two-surface, needs a composite active set) and
RoundedMohrCoulomb (not even registered — `AllYieldFunctions.h` keeps its
`#include` commented out) are OUT until P5 gives them a runtime baseline; a
specialization whose YF/PF has not opted in is refused at parse time.

**D4 — hardening is solved inside the implicit system.** `LinearHardeningFor
Scalar`, `LinearHardeningForTensor`, `NullHardening{Scalar,Tensor}` and
`ArmstrongFrederick` all become rows of `R` evaluated at `n+1`. For Null and both
Linear laws `h` does not depend on `q`, so the `q`-row is explicit and free; for
AF the `q`-row is *linear* in `alpha_{n+1}` and closes in closed form (below), so
implicitness costs one division, not an inner Newton. This is the direct removal
of M3's path dependence: BE evaluates AF's recovery term
`-c_r * ||dev m||_eq * alpha` at the *trial* `alpha` inside an otherwise implicit
loop, so its answer depends on how many Newton steps the loop happened to take.

**D5 — the explicit integrators are gated, not rewritten.** ADR-94 M2/M8 are
real (empty drift checks, unconditional `dT_min` accepts, `f_MC` in the
thousands), but rewriting four explicit integrators is a different ADR. A new
integration option `experimental_integrator 1` (default 0) is required to select
`Forward_Euler`, `Forward_Euler_Subincrement`, `Modified_Euler_Error_Control` or
`Runge_Kutta_45_Error_Control`; without it the parser refuses, cites ADR-97 and
ADR-94 M2/M8, and names the supported implicit pair (`Backward_Euler` +
`Closest_Point`). This follows ADR-94's precedent of refusing
`Backward_Euler_LineSearch` and `Runge_Kutta_45_Error_Control_old` outright, one
notch softer because these four have real users.

**D6 — `Closest_Point` evaluates `E(sigma_{n+1})` inside the residual; the
`dE/dsigma` Jacobian block is optional and defaults to zero.** Putting `E` at
`n+1` in the residual is one line and makes the converged state hyperelastically
consistent, which is the actual M9 defect: the committed stress currently
depends on which integrator evaluated `E`. Leaving `dE/dsigma` out of the
*Jacobian* is an inexact-Newton choice — it costs iterations, never accuracy,
because the converged `x` still satisfies the exact residual. That trade is free
here: **`LinearIsotropic3D_EL` is the only elasticity registered in all 43
non-StiffSoil specializations** (`DuncanChang_EL` is commented out of the
generator), and for it `dE/dsigma` is *identically zero*. The block is declared
as `ELASTICITY_STRESS_DERIVATIVE` with an inert zero base default;
`StiffSoil_EL`/`DuncanChang_EL` supply it in P5. Until they do, `tangent_type
Algorithmic` on a stress-dependent elasticity emits a one-time `opserr` warning
that `C_alg` is missing the `E,sigma : (sigma - sigma_tr)` term — that term *is*
part of the exact tangent even when dropped from the Jacobian.
`Backward_Euler` is untouched by all of this (D1).

## Residuals and Jacobians

**Conventions (fixed by ADR-94 wp/94c; do not re-derive).** Stress is
tension-positive, `p = meanStress() = trace/3`. Voigt storage is
`v = [11, 22, 33, 12, 23, 13]`; stress-like vectors store tensor shear, strain-like
vectors store *engineering* shear `gamma = 2*eps`. Consequently the plain dot of a
stress-like and a strain-like Voigt vector IS the tensor double contraction, and
`E` carries `mu` (not `2*mu`) on its shear diagonal. Every gradient in the tree
(`df_dsigma_ij`, `pf`, `df_dalpha`) is the derivative **with respect to the stored
slot**, i.e. strain-like/engineering ("Voigt convention", shear slots doubled
relative to the tensor derivative), and every contraction is a plain `.dot()`.
Write `W_s = diag(1,1,1,2,2,2)` (stress-like contraction weights, i.e.
`tensor_dot_stress_like(a,b) = a' W_s b`) and `W_e = diag(1,1,1,1/2,1/2,1/2)`
(`tensor_dot_engineering_strain_like`). New second derivatives must be produced in
the same convention: `d(n_i)/d(v_j)` with both indices on stored slots.

**Unknowns.** `x = (s, q, dl)` with `s = sigma_{n+1}` in R^6; `q = q_{n+1}` the
concatenated internal variables (a scalar per `YieldStress`/`DP_cohesion`, six
stress-like slots per `BackStress`); `dl = dlambda >= 0`, one per active surface
(1 for the smooth families and MC faces, 2 on an MC/HB edge or the MC-cutoff
corner, 3 at an MC apex or the MCTC compound corner).

**Residual.**

    R_sigma(x) = s - s_tr + dl * E(s) * m(s, q)                     (6 rows)
    R_q(x)     = q - q_n - dl * h(s, q, m(s,q))                     (n_q rows)
    R_f(x)     = f(s, q)                                            (1 row per active surface)

with `s_tr = sigma_n + E(sigma_n) * depsilon` for the predictor and
`E(s)` re-evaluated inside the residual per D6. For multiple active surfaces
`R_sigma` carries `sum_a dl_a * E * m_a` and `R_q` carries `sum_a dl_a * h_a`.

**Jacobian** `J = dR/dx`, block by block (rows: sigma, q, f; columns: s, q, dl):

    J_ss = I + dl * E * (dm/ds)            [+ dl * (dE/ds) : m , zero for LinearIsotropic3D_EL, D6]
    J_sq = dl * E * (dm/dq)
    J_sl = E * m
    J_qs = -dl * ( dh/ds + (dh/dm) * (dm/ds) )
    J_qq = I - dl * ( dh/dq + (dh/dm) * (dm/dq) )
    J_ql = -h
    J_fs = (df/ds)^T   ( = n^T )
    J_fq = (df/dq)^T
    J_fl = 0

`dh/dm` is a separate block because every hardening policy in the tree is a
function of `m`, not of `s` directly: `LinearHardeningForScalar` is
`H*sqrt(2/3 * <m,m>_e)`, `LinearHardeningForTensor` is `H*dev(m)`, AF is
`h_a*dev(m) - c_r*||dev m||_eq*alpha_dev`. Only AF has a nonzero `dh/dq`; none
has a nonzero `dh/ds`.

**Consistent tangent.** The converged system depends on `epsilon_{n+1}` only
through `s_tr`, and `d(s_tr)/d(epsilon) = E`. Differentiating `R(x(eps)) = 0`:

    J * [ dsigma ; dq ; ddl ] = [ E * deps ; 0 ; 0 ]
    C_alg = ( J^{-1} )_{sigma,sigma-block} * E              (the 6x6 top-left block of J^{-1}, times E)

Schur complement, single smooth surface. Define the **algorithmic elastic
modulus**

    Xi = ( E^{-1} + dl * dm/ds )^{-1}

(the first row of `J` rearranged: `(E^{-1} + dl*dm/ds) dsigma = deps - m*ddl -
dl*(dm/dq) dq`). Eliminating `dq` with the `q`-row and `ddl` with the
consistency row gives

    C_alg = Xi - ( Xi * m_eff ) ( n_eff^T * Xi ) / ( n_eff^T * Xi * m_eff + H_bar )

where `m_eff = m + dl * (dm/dq) * (J_qq)^{-1} * h`, `n_eff = n + (J_qq^{-T} *
df/dq)`-corrected likewise, and `H_bar = -(df/dq)^T (J_qq)^{-1} h` is the
generalized plastic modulus. With no hardening (`Null`) this collapses to the
textbook `C_alg = Xi - (Xi:m)(n:Xi)/(n:Xi:m)`; with linear hardening and
`dm/dq = 0` it collapses to `C_alg = Xi - (Xi:m)(n:Xi)/(n:Xi:m + H)` — which is
the shipped `Continuum` operator with `E` replaced by `Xi`. **That replacement
is the entire finding of M3:** setting `dl = 0` gives `Xi = E` and recovers
`Continuum` exactly, so the 57 % `Continuum` error is the size of the
`dl * dm/ds` term at Cerro-Lindo step sizes. Implementation note: form `J`
(size `6 + n_q + n_a`), LU-factor it once at convergence, and solve the six
right-hand sides `[E e_i; 0; 0]` — the Schur form above is for the *doc*; the
code reuses the one factorization it already has. `C_alg` is **unsymmetric**
whenever `m != n` (non-associated DP `etabar != eta`, MC `psi != phi`, HB
`mb_psi != mb`), so tests run on `system UmfPack`.

**Blocks that exist today vs. blocks that are new.**

| Block | Today | New member (macro name, base default) |
|---|---|---|
| `f` | `YIELD_FUNCTION` | — |
| `n = df/ds` | `YIELD_FUNCTION_STRESS_DERIVATIVE` | — |
| `m` | `PLASTIC_FLOW_DIRECTION` | — |
| `h` | `InternalVariableType::hardening_function` | — |
| `-(df/dq).h` (contracted) | `YIELD_FUNCTION_HARDENING` (`df_dxi_star_h_star`) | — (kept; CP does not use it) |
| `df/dq` (**un**contracted) | **absent** — only the contracted scalar exists | `YIELD_FUNCTION_IV_DERIVATIVE` → `df_dq(...)`, base = zero, trait `yf_has_cp_derivatives` |
| `dm/ds` (6x6) | absent | `PLASTIC_FLOW_STRESS_DERIVATIVE` → `dm_dsigma(...)`, base = one-time-warned central difference, trait `pf_has_cp_derivatives` |
| `dm/dq` (6 x n_q) | absent | `PLASTIC_FLOW_IV_DERIVATIVE` → `dm_dq(...)`, base = zero |
| `dh/dq`, `dh/dm` | absent | `HARDENING_FUNCTION_DERIVATIVES` → `df_dq`, `df_dm` statics on each policy, base = zero |
| `dE/ds` | absent | `ELASTICITY_STRESS_DERIVATIVE`, base = zero (D6) |

`dh/ds` is not declared: no policy in the tree reads `sigma` (they read `m`,
`depsilon` and `current_value`); a policy that later does must declare it.

**VonMises (analytic).** `f = ||r||_s - sqrt(2/3)*k`, `r = dev(s) - alpha`,
`||r||_s = sqrt(r' W_s r)`, `n = m = W_e^{-1} (I_dev' W_s r)/||r||_s` in the
stored/Voigt convention (this is exactly the shipped `df_dsigma_ij` with its
three doubled shear slots). Then

    dm/ds  = ( I_dev - nhat (x) nhat ) / ||r||_s        (Voigt-weighted; nhat = r/||r||_s)
    dm/dalpha = -dm/ds  (restricted to the deviatoric subspace)
    df/dk  = -sqrt(2/3)        df/dalpha = -n
    h_k    = H * sqrt( (2/3) <m,m>_e )        dh_k/dm = (2H/3) * W_e m / sqrt((2/3)<m,m>_e)
    h_alpha (Linear) = H_kin * dev(m)         dh_alpha/dm = H_kin * I_dev,   dh_alpha/dalpha = 0

Associated flow, so `dm/ds` is the yield Hessian and `C_alg` is symmetric.
Radial-return closed form is available and is the P1 oracle
(`adr94_oracle/vm_shear_oracle.py` already pins the shipped map at 1.5e-13).

**Drucker–Prager (analytic).** `f = sqrt(J2) + eta*p - xi_c`,
`sqrt(J2) = sqrt(0.5 * r' W_s r)`, `g = sqrt(J2) + etabar*p`. Both gradients are
the shipped ones (normal slots halved, `eta/3` on the pressure part). Hessian:

    d(sqrt J2)/ds ds = ( I_dev / (2 sqrt J2) ) - ( r (x) r ) / ( 4 * J2^{3/2} )     (Voigt-weighted)
    dm/ds = the same expression with etabar (the eta*p term has zero Hessian)
    dm/dalpha = -dm/ds|_dev ;  df/dalpha = -(1/2)-weighted-normal-slot form (shipped df_dalpha)
    df/dxi_c = -1 ;  h_xi_c = H * sqrt(2/3 <m,m>_e)   (LinearHardeningForScalar)

**Drucker–Prager apex.** The vertex is a measure-zero point of `f` where
`dm/ds` blows up as `1/sqrt(J2)`, so it gets its own reduced system. With
`dev(sigma_{n+1}) = 0`, `x` collapses to `(p, dl)`:

    R_p = p - p_tr + dl * K * etabar_vol      R_f = eta*p - xi_c(q) = 0
    => p_{n+1} = xi_c/eta ,  d(eps^p)_vol = ( p_tr - xi_c/eta ) / K

and the consistent tangent is the **bulk-only rank-1 projector**

    C_alg = beta * (1 (x) 1) / 3 ,  1 = [1,1,1,0,0,0]',  beta = 0 for perfect plasticity,
    beta = d(xi_c)/d(eps^p_vol) / eta  when xi_c hardens.

The region test is where CP improves on the shipped code: `check_apex_region`
is a **Euclidean** normal-cone test `p - p_apex >= eta * q` (the honest caveat is
in its own comment and in `LEDGER_quirks`), while the exact condition is in the
**elastic metric**, `p_tr - p_apex >= (K*etabar/G) * q_tr`. `Backward_Euler`
cannot fix that (the YF signature cannot see `K`, `G`); `Closest_Point` can,
because classification happens in the integrator where `E` is in scope. **CP does
its own region classification and does not call `check_apex_region`.** The
`f(sigma_apex) <= tol` admissibility gate is kept, as wp/94c wrote it.

**MohrCoulomb (principal space).** Work in tension-positive principal stresses
`(sig1 >= sig2 >= sig3)` from `sigma.principalStresses()` (which returns
ascending). In that space MC is six planes; on the sorted sextant it is the
single plane `f = (sig1 - sig3) + (sig1 + sig3) sin(phi) - 2 c cos(phi)`, and the
flow potential the same with `psi`. The return is then a **linear** projection in
the elastic metric onto (a) the plane, (b) one of two intersection lines
(triaxial-compression and triaxial-extension edges), or (c) the apex point
`p_apex * 1`. Region selection is Clausen's **boundary-plane test on the trial
state**: each boundary plane is spanned by the intersection line's direction and
the flow direction of the adjacent surface, and the trial principal point's side
of it decides the region — an exact, branch-free classification that replaces the
`|theta| < 29 deg` guard entirely. For an edge or the apex, 2 or 3 surfaces are
active and the plastic multipliers follow Koiter's rule; the principal-space
consistent tangent is the closed form of Clausen, Damkilde & Andersen (2007),
constant on each region.

Transforming back to 6D is where the work is. With `A` the matrix of eigen-
projections of the **trial** stress,

    C_alg(6x6) = A^T * C_princ * A  +  sum_{i<j} [ (sig_i^ret - sig_j^ret) / (sig_i^tr - sig_j^tr) ] * S_ij

where `S_ij` are the spin/rotation contributions from the derivative of the
eigenprojections. The bracketed ratio is the standard 0/0 at a **degenerate
trial eigenvalue** (`sig_i^tr == sig_j^tr`, i.e. any axisymmetric path — which
is most triaxial decks); the limit is `d(sig_i^ret)/d(sig_j^tr)` taken from
`C_princ`, so the implementation switches to the l'Hôpital value when
`|sig_i^tr - sig_j^tr| < eps_deg * strength_scale()`. That relative threshold is
not optional: an absolute one reintroduces M5's unit dependence.

**MohrCoulombTensionCutoff.** One more plane, `f_TC = sig_max - T`, i.e. one more
active surface; the compound `MC ∩ TC` corner is a 3-surface Koiter return.
ADR-84's `SPECIAL_RETURN` already computes exactly this geometry (cutoff face,
Rankine edge, `MC ∩ TC` corner, apex `T_eff * 1` with
`T_eff = min(TC_min_stress, c*cot(phi))`) and already returns the **raw
active-set (Koiter) tangent** in `stiffness_return` with an `SR_QUALITY_*` flag
— P3 of that ADR stopped blending inside the hook precisely so the integrator
could apply `tangent_type`. `Closest_Point` reuses the hook verbatim, maps
`SR_QUALITY_EXACT` to `Algorithmic = stiffness_return` (it *is* the consistent
tangent of that exact return), and treats `SR_QUALITY_FALLBACK` as a refusal
under `strict_convergence`, as `Backward_Euler` does. No ADR-84 geometry is
re-derived here.

**HoekBrown.** `f = max(f_shear, f_tension)` (jaabell's composite, ported in
#806) in compression-positive principal stresses; `f_shear = sig1 - sig3 -
sigma_ci*(mb*sig3/sigma_ci + s)^a`, `f_tension = sigma_t - sig3`. Both shipped
gradients (`HoekBrown_YF::df_dsigma_ij`, `HoekBrown_PF`) are **central
differences over the six raw Voigt slots** — the header says so: "the analytical
derivative involves principal stress directions which can be numerically
sensitive". In principal space that objection evaporates: `df/dsig1 = 1`,
`df/dsig3 = -1 - a*mb*(mb*sig3/sigma_ci + s)^(a-1)`, `df/dsig2 = 0`, and the
same with `mb_psi` for `g`. So **HoekBrown gets analytic principal-space
derivatives, not a finite-difference fallback**; the 6D FD is retained only for
`Backward_Euler`, unchanged. The only genuine FD user is the base
`PLASTIC_FLOW_STRESS_DERIVATIVE` default, which emits a one-time
`opserr` warning naming the PF class; no family shipped by this ADR should ever
reach it. Note also a pre-existing mismatch to fix here: `HoekBrown_PF::g` uses
an `if (arg > 0) else` branch while the YF uses the composite `max` — the two
disagree in the tension regime.

**P0 measured all of this** (`adr97_oracle/cppm_hb.py`, README HB block); four
results change what P3 has to build. (i) The `g`/`f` mismatch is *not* confined
to tension: `HoekBrown_PF::g` never negates to the compression frame, so its
branch `arg` is built from the most COMPRESSIVE principal, is negative on every
compressive state, and `g` collapses to a **Tresca** potential — non-dilatant,
independent of `HB_mb_psi`, 32.86° off the normal even at `mb_psi = mb`. At the
apex all six shipped flow directions have negative trace, so a trial past the
tensile corner has **no** return to the apex at all — the mechanism behind the
ADR-94 H10a residual. P3 fixes the frame; it is not a tension-only edit.
(ii) The composite's tension branch is **inert on the yield surface** (`f_shear
≤ 0` already forces `sig3 >= sigma_t`), so `Closest_Point` needs **no
tension-plane return** — only the apex, which is the shear/tension corner.
(iii) `CHECK_APEX_REGION` is Euclidean where the exact region is
`apex + D3·(positive octant)`; since `D3·octant` is a strict subset it always
OVER-claims (32/400 trials, up to 122.6 kPa of silently lost strength).
(iv) Two formulation requirements: the Newton must use the surface's natural
variable (`arg = w^(2/a)`; with `sig3`/`y1` as the unknown it leaves the domain
on step 1 and the next Jacobian is singular) and must **normalize the flow
direction** (`|m| ~ arg^(a-1)` blows up at the apex) — together these are what
hold the `<= 5` iteration gate. Near the apex `|f|`'s own round-off floor
exceeds `1e-10`, so the yield tolerance must be gradient-scaled.

**ArmstrongFrederick, implicit (this is what removes M3).**

    h_alpha = h_a * dev(m) - c_r * ||dev m||_eq * alpha_dev ,   ||v||_eq = sqrt( (2/3) <v,v>_e or _s )
    alpha_{n+1,dev} = ( alpha_{n,dev} + dl * h_a * dev(m) ) / ( 1 + dl * c_r * ||dev m||_eq )

The `q`-row is **linear in `alpha_{n+1}`**, so it closes in closed form and the
implicit update costs one division. Blocks:

    dh_alpha/dalpha = -c_r * ||dev m||_eq * I_dev
    dh_alpha/dm     = h_a * I_dev  -  ( (2/3) * W_e * dev(m) / ||dev m||_eq ) (x) alpha_dev
    J_qq            = I + dl * c_r * ||dev m||_eq * I_dev        (well-conditioned, always)

Two consequences. (1) The closed form is a **contraction toward the saturation
sphere** `||alpha|| = h_a/c_r`, so `alpha_{n+1}` cannot overshoot it — the
shipped hard `if (alpha_norm >= alpha_limit) derivative = 0` branch is both
unnecessary inside CP and non-differentiable, so `Closest_Point` drops it
(`Backward_Euler` keeps it, D1). (2) AF's `h` reads `depsilon`, the step
increment; CP passes the same argument, so the signature is unchanged.

**The Newton loop.**

- Start from the elastic predictor: `s = s_tr`, `q = q_n`, `dl = 0`.
- Convergence: `|R_f| <= tol_f` with `tol_f = max(f_absolute_tol,
  f_relative_tol * yf.strength_scale())` — the existing `yf_tolerance()`
  accessor, unchanged (M5) — **and** `||R_sigma|| <= max(stress_absolute_tol,
  f_relative_tol * strength_scale())` **and** a matching scaled norm on `R_q`
  (each IV normalized by its own `strength_scale`-derived reference, so a
  `BackStress` in Pa and a `YieldStress` in kPa do not fight).
- `n_max_iterations` bounds the loop; a converged CPPM from an elastic predictor
  is quadratic and should take <= 5.
- **No line search for the smooth families.** For convex `f` with associated
  flow the CPPM residual is the stationarity system of a strictly convex
  closest-point projection in the `E`-metric, so Newton from the elastic
  predictor is globally convergent on the yielding branch; non-associated, the
  map is still a contraction in a neighbourhood whose radius scales with
  `||E^{-1}||/||dm/ds||`, and implicit-FE step sizes sit inside it. A line search
  that "succeeds" on a non-converged system is ADR-94 M7's
  `Backward_Euler_LineSearch` failure mode (2/20 vs BE's 20/20).
- Failure: `LADRUNO_MATERIAL_REFUSED` on non-convergence, on NaN in `x` or `R`,
  and on a singular/ill-conditioned `J` (LU failure, or reciprocal condition
  below a fixed floor). Never a bare `-1` (ADR-94 B2: `LadrunoBrick` compares
  only against the sentinel; `stdBrick` drops everything).
- Every commit path ends in the existing `ladruno_strict_rejects(...)` gate, so
  `strict_convergence 1` refuses `|f| > tol` from `Closest_Point` on the same
  contract it already has for `Backward_Euler`.

**`Closest_Point` does NOT reuse `Backward_Euler`'s intersection/elastic-fraction
split — and neither does BE.** Confirmed by reading: `Backward_Euler` calls
`intersection_stress.setZero()` / `intersection_strain.setZero()` at setup and
**never reads either again**; the yield-crossing bisection lives in
`Forward_Euler`, `Forward_Euler_Subincrement` and `compute_local_stress()` only.
So the split affects BE's committed state not at all, and CP dropping it changes
nothing relative to BE. What *does* differ is the map itself (cutting plane vs
closest point). Consequence for gate 4, stated honestly: **CP and BE agree to
Newton tolerance exactly when the flow direction does not rotate over the step**
(von Mises proportional loading, DP proportional loading, any radial path), because
then `sum_k dl_k m(sigma^k) == dl * m(sigma_{n+1})`. On rotating-normal and on AF
decks they differ by the cutting-plane path error, which is gate 3's pinned
contrast, not a gate-4 failure. The gate-4 "perfectly plastic decks" must
therefore be non-rotating-normal decks, named as such in the test.

## Phases

Every WP: its own `wp/97<letter>-<slug>` branch cut from fresh `ladruno`, its own
worktree, a **day-one draft PR**, commits early and often, **one build per WP**
reusing the warm build cache in `../asdplastic-review-plan-62585c` (do NOT create
a fresh MUMPS build), and every source edit staged as a **re-runnable script** so
a session restart does not orphan the work (ADR-94 lesson; WMI-launched builds
survive, `Start-Process` ones do not).

| P | WP branch | Content | Model | Effort |
|---|---|---|---|---|
| **P0** | `wp/97a-plan-oracles` | This plan + the numpy/g++ oracles in `Ladruno_implementation/adr97_oracle/`: closed-form radial return (VM), DP flank + apex, MC principal-space face/edge/apex with the Clausen tangent, AF one-step implicit update, and a 6x6 `C_alg` reference from a dense FD of each oracle's own committed map. No C++. | opus | 1 session |
| **P1** | `wp/97b-cp-smooth` | `Closest_Point` + `Algorithmic` for VonMises and Drucker–Prager (flank + apex), Null/Linear-scalar/Linear-tensor/AF hardening; the six new interface macros with inert base defaults; parser tokens `Closest_Point`, `Algorithmic` + the D2 cross-refusal. **20 of the 46 specializations** supported (IV_YF: VM 2 + DP 2; IV_PF: VM_PF 3 + DP_PF 2). | opus | 2 sessions |
| **P2** | `wp/97c-cp-principal` | Principal-space multi-surface CPPM + Koiter tangent for MohrCoulomb, and MohrCoulombTensionCutoff on top of ADR-84's `special_return`; the eigenprojection back-transform incl. the degenerate-eigenvalue limit. Support rises to **31** (YF {VM,DP,MC} x PF {VM,DP,MC} = 30, + MCTC). | opus | 2 sessions |
| **P3** | `wp/97d-cp-hoekbrown` | HoekBrown principal-space return to the CURVED surface (Clausen & Damkilde 2008) with analytic principal derivatives; the `HoekBrown_PF::g` / composite-`f` tension mismatch. **Support rises to 23 of 46, NOT the 43 written here.** Same correction as P2's: the matched-pair rule is what makes the map verifiable, and `HoekBrown_YF` x `HoekBrown_PF` is the ONLY registered pairing with both functors of the family (`HoekBrown_YF` appears in 7 specializations, `HoekBrown_PF` in 6). The other eleven stay refused. | opus | 1 session |
| **P4** | `wp/97e-numalg-repoint` | Re-point `Numerical_Algorithmic_FirstOrder/SecondOrder` at the actual committed map (closes ADR-84 P4 and M3's "third map"); `compute_local_stress()` retained only as a private helper, marked as not-a-map. | sonnet | 1 session |
| **P5** | `wp/97f-explicit-gate` | D5's `experimental_integrator` gate + refusals; first runtime coverage for `StiffSoilShear`/`StiffSoilCap` (M9) and the root cause of the `StiffSoilShear` step-1 NaN (its PF is non-finite at 6/194 cloud points); `ELASTICITY_STRESS_DERIVATIVE` for `StiffSoil_EL`/`DuncanChang_EL` if the coverage lands. | sonnet | 1–2 sessions |
| **P6** | `wp/97g-measure` | Cerro-Lindo-scale measurement: wall-clock and `testIter()` for {BE,CP} x {Secant,Continuum,Algorithmic,Numerical_Algorithmic} on a real mesh; the input to the default-flip decision (a separate PR, D1). | sonnet | 1 session |
| **P7** | `wp/97h-closeout` | Verdict + verification manifest + mutation-gate record + user guide (ADR-87 warrant package); banner line; `LEDGER_vanilla_files` rows, `LEDGER_implementations` row, `LEDGER_quirks` entries; this doc's implementation log. | sonnet | 1 session |

The remaining 23 -- the eleven mixed HoekBrown pairings, the six mixed MohrCoulomb
pairings, the three VonMises/DruckerPrager x MohrCoulomb crosses counted among them,
and the StiffSoil trio -- are **refused at parse time** under
`Closest_Point` unless P5 delivers; the refusal names the family and cites D3.

## Gates

1. **Material-level correctness, per family.** `|f| <= tol` at *every* commit;
   committed stress vs oracle `<= 1e-10` relative on triaxial, simple-shear and
   rotating-normal paths; the local Newton converges quadratically in `<= 5`
   iterations.
2. **Tangent.** `Algorithmic` vs a central difference of the material's OWN
   committed response `<= 1e-6` relative on a FREE-DOF rig, and
   `ops.testIter() <= 3` per step on the ADR-94 two-cube heterogeneous model
   (`tests/test_adr94_redblue_blue.py`) against the recorded `Continuum` count.
3. **Path independence.** Step-halving reproduces the AF / linear-hardening state
   to `<= 1e-9`; the cutting plane fails this by construction — pin the contrast.
4. **Default inertness.** `Backward_Euler` + `Secant` histories byte-identical
   before/after, in FRESH SUBPROCESSES; CP vs BE agree on perfectly plastic decks;
   CP vs BE differ measurably on AF decks.
5. **Mutation gate.** Drop one consistent-tangent term (e.g. the `dl * dm/ds`
   term of `Xi`) on a scratch build; the gate-2 tests go red; recorded in the
   manifest.
6. **Fail-loud.** `LADRUNO_MATERIAL_REFUSED` on every new failure path; strict
   mode refuses `|f| > tol`; the parser rejects unknown tokens.
7. **Portability.** Zone-A green on Linux; cross-platform float pins
   `>= 1e-6` relative.

**Mapping.**

| Gate | Test | Oracle |
|---|---|---|
| 1 (VM, DP, apex) | `tests/test_adr97_p1_smooth.py` | `adr97_oracle/cppm_vm.py`, `cppm_dp.py` |
| 1 (MC, MCTC) | `tests/test_adr97_p2_principal.py` | `adr97_oracle/cppm_mc.py`, ADR-84 `test_asdplastic_mctc.py` rows |
| 1 (HB) | `tests/test_adr97_p3_hoekbrown.py` | `adr97_oracle/cppm_hb.py` (**written in P0**; see the README's HB block for the pins and for the four new header findings it measures) |
| 2 | `tests/test_adr97_p1_smooth.py::test_algorithmic_vs_fd`, `..._p2_principal.py` (same name) | `adr97_oracle/fd_tangent_driver.py` (`fd_check()`; BE/Continuum negative control 0.573) + the recorded `Continuum` iteration count from `test_adr94_redblue_blue.py` |
| 3 | `tests/test_adr97_p1_smooth.py::test_af_step_halving` | `adr97_oracle/path_independence.py` + `cppm_vm.py` case (c) |
| 4 | `tests/test_adr97_p4_inertness.py` (fresh-subprocess helper `_run_child`, ADR-94 pattern) | pre-change stress dumps, same binary |
| 5 | `tests/test_adr97_p5_mutation.py` + the manifest's mutation record | — (scratch build) |
| 6 | `tests/test_adr97_p6_faillloud.py` | — (stderr assertions; tet or `LadrunoBrick` host only) |
| 7 | Zone-A CI (`ladruno.yml`) | — |

Every gate-1/2/3 test runs on `TenNodeTetrahedron` or `LadrunoBrick` with
`system UmfPack`. `stdBrick` is used only as a *negative* control: it swallows
every material return code by design (ADR-94 B2, pinned).

## Risks / open questions

> [!question]
> Does the principal-space MC return stay quadratic at the sextant boundary
> `theta = +-30 deg`, where the return region changes discontinuously? Clausen's
> boundary-plane test is exact, so the *classification* does not chatter — but a
> trial state that crosses a boundary plane between global Newton iterations can
> make the element-level residual non-smooth. Measure on gate 2's iteration count.

> [!question]
> Is `dE/dsigma` (D6) worth implementing for `StiffSoil_EL` at all, or should
> stress-dependent elasticity simply refuse `tangent_type Algorithmic`? Decide at
> P5 with data, not now.

- **MSVC != GCC on temporaries.** Never bind an Eigen product temporary to a
  non-const reference (`const VoigtVector& Em = Eelastic * m;`) — compiles on
  MSVC, dangling on GCC; it cost the ADR-94 wave a red CI run. Named locals only.
- **`/bigobj` TU rebuild cost.** `OPS_AllASDPlasticMaterial3Ds.cpp` instantiates
  all 46 specializations in one TU; any header edit rebuilds it, and six new
  members per YF/PF grows it further. Use
  `ninja -j1 CMakeFiles/OPS_Material.dir/.../OPS_AllASDPlasticMaterial3Ds.cpp.obj`
  for syntax turnaround (`LEDGER_quirks`).
- **Per-tag static option maps.** `INT_OPT_*`/`DBL_OPT_*` are
  `std::map<int, ...>` keyed by `ASDP_TAG` — shared across every instance of a
  tag *by design*. The new tokens follow that pattern; do not "fix" it, and do
  not add a per-GP `map::operator[]` lookup inside the Newton loop (ADR-94 m1).
- **Euclidean vs elastic-metric apex classification** (`LEDGER_quirks`, wp/94c):
  CP fixes it for DP by classifying inside the integrator, but the *shipped*
  `check_apex_region` stays Euclidean for `Backward_Euler`. Two different answers
  for the same YF on two integrators is a documentation obligation, not a bug.
- **Unsymmetric tangents.** `C_alg` is unsymmetric for every non-associated
  family. Tests use `system UmfPack`; `ProfileSPD` is wrong and `FullGeneral`
  crashes fully-prescribed material-point rigs (`FullGenLinSOE` N=0,
  `LEDGER_quirks`). PARDISO's symmetric `-matrixType` (ADR-75 P1d) must not be
  selected for these models.
- **Eigen-decomposition cost per Gauss point.** MC/MCTC/HB now do a 3x3
  symmetric eigensolve (with eigenvectors) per return, where the 6D map did an FD
  gradient loop of 12 `f` evaluations. Likely a wash or a win — a P6 measurement,
  not an assumption.
- **Degenerate eigenvalues.** The `(sig_i^ret - sig_j^ret)/(sig_i^tr - sig_j^tr)`
  spin term is 0/0 on every axisymmetric path — most triaxial decks, so this is
  the common case, not the corner case. The l'Hôpital branch needs its own test
  and a threshold relative to `strength_scale()` (M5).
- **Backward compatibility.** D1 keeps `Backward_Euler` byte-identical, but gate
  4 must run in **fresh subprocesses**: the per-tag option maps and the pytest
  heap both leak state across tests in one process (ADR-94 `capfd`/`WinError 6`).
- **`Algorithmic` on the apex/`special_return` paths.** Both in-`Backward_Euler`
  tangent switches already list `TOT::Algorithmic` alongside `Continuum`; those
  cases become unreachable for BE once D2's refusal lands. Leave them, marked.

## Implementation log

- **2026-09-07** — plan drafted (opus), draft PR #817 on `wp/97a-plan-oracles`. P0 oracles
  delivered the same day (opus, `adr97_oracle/`, `reference_output.txt` on build `3622d6214`):
  FD checks 9.6e-12 … 1.7e-10 on every consistent tangent (complex-step Jacobians); the
  element-level free-DOF driver reproduces ADR-94 M3 from the binary alone (BE/Continuum
  0.573, Secant 0.799, Elastic 1.025, Numerical_* 0.046); CPPM per-step invariance to the
  Newton start 1.1e-12 vs a cutting-plane spread of 9.15 (σ) / 5.68 (α) on a ~46 stress for
  VM+AF — with LINEAR hardening the spread is 7e-15, which is why ADR-94 H6 could not see it.
- **P0 header findings that P1 must decide on (see `adr97_oracle/README.md` §"Header findings"):**
  (1) `ArmstrongFrederickPolicy` adds the engineering-shear `mdev` to the stress-like back
  stress (shear slots grow 2×) and has no 2/3 on `ha` — the oracle pins the HEADER convention;
  the policy is shared with `Backward_Euler`, so under D1 `Closest_Point` mirrors it and the
  convention question is recorded, not fixed. (2) `DruckerPrager_YF` has its cohesion IV
  commented out of `f` while `yf_hardening` still contributes `df/dk = −1`: `Closest_Point`'s
  uncontracted `df/dq` must be the TRUE derivative of the header's `f` (zero for that IV), so
  DP-with-cohesion-hardening is perfectly plastic under CP and hardening under BE until the
  YF is fixed (separate PR; it changes BE). (3) DP/MC apex-region tests are Euclidean, not
  elastic-metric — quantified misclassifications in the oracle; CP must use the elastic-metric
  test (its own code path, BE untouched). (4) `MohrCoulomb_PF` converts `MC_c` to radians
  (inert). (5) MC has no edge/apex algebra; `cppm_mc.py` is the reference. (6) At the MC apex
  with ψ<φ the Koiter multipliers are not all positive: choose the region by boundary planes,
  never by a dΛ≥0 active-set search.

- **2026-09-07 — P1 (`wp/97b-cp-smooth`, PR [#819](https://github.com/nmorabowen/OpenSees/pull/819), build `dd05b60a6` then `ec6091c4f`).** `integration_method Closest_Point` + `tangent_type Algorithmic` SHIPPED for VonMises and Drucker–Prager (flank + apex) with Null / LinearScalar / LinearTensor / ArmstrongFrederick hardening — 20 of 46 specializations. Six new interface members with inert base defaults + opt-in traits; the parser refuses `Closest_Point` for unconverted families and `Algorithmic` for every other integrator (D2). **Measured:** gate 1 — committed stress vs the P0 oracles 7.4e-14 … 2.4e-12 (VM perfect/linear, 3 paths), 1.3e-13 / 9.9e-13 / **9.2e-9** (VM+AF triaxial / simple shear / rotating normal; the 9.2e-9 is the yield tolerance accumulated over 20 steps, `f_absolute_tol` 1e-6 on a ~45 stress), DP cone associated 4.9e-13 and non-associated 5.1e-13, DP apex exact; worst committed `|f|` over a whole path 7.1e-15; local Newton **1** iteration (perfect, linear) and **3** (AF), gate `<= 5`. gate 2 — `Algorithmic` vs a central difference of the binary's own assembled residual **1.51e-11** (4 free DOFs) and **2.15e-10** (12-DOF sheared rig) against the pinned `Backward_Euler`/`Continuum` **0.573447** reproduced to the digit; ADR-94 two-cube `testIter` **4,4,4,4** vs `Continuum`'s **6,41,36,30** (7.1× fewer; the ADR's `<= 3` was an estimate, 4 is what an exact tangent costs on that rig). gate 3 — AF step-refinement error at N=10 against each map's own N=160 limit: CP 5.68e-3 vs BE 2.31e-2, **4.07×** (P0 predicted 4.3×); with LINEAR hardening the two maps agree to **1.45e-16**, the reason ADR-94 H6 was blind to this. gate 4 — 23 decks / 282 rows byte-identical in fresh subprocesses; CP≡BE at 0.0 / 1.6e-33 / 3.8e-17 on non-rotating perfectly plastic decks and 2.3e-2 apart on AF. gate 6 — 13/13 loud. **Two P0 findings honoured, not fixed** (both would change `Backward_Euler`, which D1 forbids): DP's cohesion IV is out of its own `f`, so DP-with-cohesion-hardening is perfectly plastic under CP and hardening under BE (pinned as a test); and AF's saturation branch is KEPT, with CP differentiating the branch `f` actually takes, because `f` is shared with BE. **Found while testing:** the elastic-metric apex test is a SIGN test on `dot(dev_ret, dev_tr)` and degenerates to "CONE" for a trial state exactly on the hydrostatic axis — ADR-94 B4's own reproducer; fixed with a `||dev_tr|| <= tol_f` short circuit. Also `dm_dsigma_buffer` needed `this->` (dependent base, GCC-only error) and a `Path` time series returns 0 at its last time point. **gate 5 (mutation, [[reviews/adr97_p1_mutation]]):** dropping the `dl*dm/ds` term of `Xi` from the Jacobian on a scratch build turns `Algorithmic` into EXACTLY `Continuum` -- the free-DOF FD error goes 1.51e-11 -> **0.573447** (4 DOF) and 2.15e-10 -> **0.670133** (12 DOF), both matching the P0 oracle's `Continuum` values to every digit, plus the AF local Newton 3 -> 13 iterations; 4 killed, 39 survivors, each one a case the mutation provably cannot reach (the residual is untouched, so the committed stress stays exact). Reverted and re-measured identical. **The two P0 findings above also SHARPENED one of them:** measured, the DP cohesion hardening is inert on the committed stress under BOTH integrators (CP exactly 0, BE 4.5e-9 of cutting-plane iterate-path residue) -- it cannot be otherwise, since both maps drive the SAME `f` to zero and `f` does not read the IV; what differs is the TANGENT, because BE's `Continuum`/`Secant` divide by `n:E:m - H` with that phantom `H`. The IV itself still grows (k = 0.190036) and is read by nothing.
- **2026-09-07 — P2 (`wp/97c-cp-principal`, PR [#824](https://github.com/nmorabowen/OpenSees/pull/824), build `1072c27ae`).** `integration_method Closest_Point` + `tangent_type Algorithmic` SHIPPED for the **Mohr-Coulomb family** through a PRINCIPAL-STRESS-SPACE multi-surface closest-point return (Clausen, Damkilde & Andersen 2006/2007): on the sorted sextant the surface is one PLANE, two corner LINES and a vertex, so every return is a closed-form linear projection in the elastic metric — no Newton at all, and `cp_iterations` reads 1 on a plastic step and 0 on an elastic one. Region selection is Clausen's BOUNDARY-PLANE test, whose signs are taken analytically (`sgn_L = sign(n_L·ℓ_other)`) rather than from the oracle's dimensional reference point — verified identical over 60 (φ,ψ,c) combinations and 13094 trial states, **0 mismatches**; the edge directions come from a cross product instead of an SVD (2.2e-16) and `Rs⁻¹` is built as the Voigt image of `Qᵀ·Q` instead of a numerical inverse (2.4e-15). **MohrCoulombTensionCutoff reuses ADR-84's `special_return` verbatim** and only falls back to the plain-MC return when that hook declines, re-checking the COMPOSITE `f`. **Measured:** gate 1a — all six P0 oracle trial states (face, sheared face, both corner lines, two apex states) to **1.5e-16 / 9.1e-16 / 0.0 / 8.0e-16 / 0.0 / 0.0** relative with committed `|f| <= 5.9e-14`; gate 1b — worst committed `f` **1.2e-14 / 1.8e-15 / 2.0e-14** on triaxial, simple-shear and rotating-principal-direction paths (the rotation checked at **20.70°**); gate 1c — MCTC hydrostatic tension `24.7·I` and the Rankine face return both **BIT-IDENTICAL to `Backward_Euler` (gap 0.0)**, which is the check that ADR-84's geometry is reused rather than re-derived, and the confined-compression fall-through admissible on BOTH branches (`f_MC` 2.8e-14 on a scale of 94, `f_TC` −217); gate 2 — `Algorithmic` vs a central difference of the binary's own assembled internal force **2.88e-11** (degenerate edge, l'Hôpital branch) and **8.40e-9 / 1.08e-8** (face, separations 0.19 / 0.27) on a NEW free-node rig, with the shipped `Backward_Euler`/`Continuum` and `/Secant` failing to converge on that same rig (`analyze -> -3`), and global-Newton **16 vs 62 / 50** (MC oedometric) and **22 vs 74** (ADR-84 MCTC); gate 4 — the 23-deck / 282-row `Backward_Euler` baseline still byte-identical, CP step-size independent (N = 1/4/10/40 all 1e-16) and CP ≡ BE at the apex to **2.1e-16**; gate 6 — 12 refusal tests, including all six MIXED YF/PF pairings, `MC_phi == 0`, and three hydrostatic-plus-vanishing-deviator trials that commit a finite vertex instead of NaN. **Support 20 → 22 of 46, not the plan's 31:** `MohrCoulomb_YF` is registered in 7 specializations and `MohrCoulomb_PF` in 6, and in only ONE are both of the family; the principal map assumes BOTH the surface and the potential are piecewise linear and the smooth 6D map cannot use MC's Lode-angle gradient, so the six mixed pairings stay refused, enforced by matching compile-time family markers plus an all-IVs-inert fold. **Found while implementing:** (1) `Backward_Euler` reproduces the exact return ONLY through its own finite difference — `MC_ds = 1e-4` gives **2.1e-14**, `MC_ds = 0` (the shipped ANALYTIC Lode-angle branch, which every deck in this repo uses) gives **2.9e-1** and is step-size dependent; pinned in both directions, NOT fixed, because it changes `Backward_Euler` (D1); (2) an oedometric MC deck at `nu = 0.25`, `phi = 30` NEVER yields, because `K0 = nu/(1-nu) = 1/3` coincides exactly with the compression meridian — a silently vacuous tangent gate; (3) `yf_tolerance()` is not a usable admissibility tolerance for a stress reassembled from a spectral decomposition (round-off 6e-10 relative refused a valid MCTC step, amplified because an edge return lands on a Lode-angle corner) — replaced by an EXACT principal-space check plus a stress-relative composite-`f` guard; (4) neither `fd_tangent_driver` rig can reach a Mohr-Coulomb FACE state. All four are in [[LEDGER_quirks]]. **gate 5 (mutation, [[reviews/adr97_p2_mutation]]).** Full report: [[reviews/adr97_p2_report]].
- **2026-09-07 — P3 (`wp/97d-cp-hoekbrown`, PR [#825](https://github.com/nmorabowen/OpenSees/pull/825), build `d0de2d7c8`).** `integration_method Closest_Point` + `tangent_type Algorithmic` SHIPPED for the **Hoek-Brown family** as a principal-stress-space return to a **CURVED** surface (Clausen & Damkilde 2008). P2's principal machinery — spectral decomposition, eigenprojection back-transform with the rotation term and its l'Hôpital limit, analytic `Rs⁻¹`, tangent policy, plastic-strain convention, strict-convergence contract — is REUSED verbatim; only the projection is replaced, by a 4×4 Newton on the face or either curved edge plus a closed-form vertex. Three things the P0 oracle measured as REQUIREMENTS: (a) the Newton runs in the surface's own variable `arg = (w²)^(1/a)` — with `y1` as the unknown the first step overshoots to `arg = −2.2e-03` and the next Jacobian is singular; (b) the flow direction is NORMALIZED in the residual, because `|m| ~ arg^(a−1)` blows up exactly where the near-apex returns land (460.6 vs 3.9), which is what holds the worst case at **5** Newton iterations instead of 6; (c) the admissibility tolerance is GRADIENT-scaled, because `|df/dy1| = 1 + a·mb·arg^(a−1)` diverges at the apex and an absolute 1e-10 gate is unattainable within ~1e-2 kPa of the vertex. Region selection uses NO boundary planes (a curved surface has none): the apex by an exact elastic-metric cone test in dual/facet form, then the face, then the edges ordered by the FACE return's own `y1−y2` / `y2−y3` margins, which ARE the exact signed boundary functions. The composite's tension branch is INERT on the surface (`f_shear ≤ 0` already forces `y1 ≤ T`; 0 of 4000 sampled surface points have `f_tension` winning, `min(f_shear − f_tension) = 1.004e-04`), so the shear/tension corner IS the apex and no Rankine face return exists — a Rankine return from the oracle's `[645,265,255]` trial leaves `f_shear = +123.34` kPa, inadmissible. **Measured:** gate 1a — all SEVEN oracle regions to **1.448e-15 / 1.735e-15 / 7.259e-16 / 1.268e-15 / 0.0 / 0.0 / 1.851e-13** relative with `cp_iterations` **4/4/3/3/1/1/5** and committed `|f| ≤ 1.8e-11`; the last row is the trial the header's EUCLIDEAN `CHECK_APEX_REGION` calls an apex, where `APEX_STRESS` would discard **122.6 kPa = 2.29 %** of the strength scale (the header always over-claims: `D3·(octant)` is a strict subset of the octant, and 32 of the oracle's 400 scanned trials disagree). gate 1b — worst committed `f` **0.0 / +7.3e-12 / +7.3e-11** on triaxial / simple-shear / rotating paths (rotation checked at **16.80°**). gate 1c — the ADR-94 uniaxial tension plateau lands on the **exact closed-form limit 244.4854419 kPa (1.16e-16 relative)**, NOT on the apex `T` = 245.0152 that ADR-94's `rel = 5e-2` pin accepts (`T` is 0.216 % above the limit, because at `y1 = T` the clamp leaves `f_shear = y1 > 0`); and a trial 2×/20× past the tensile corner now commits the finite apex instead of stalling, because the frame-consistent potential HAS a return to the vertex where the shipped Tresca `g` has none. gate 2 — `Algorithmic` vs a central difference of the binary's own assembled internal force **3.342e-09 / 3.869e-09** (free-node FACE rig, curvature AND rotation terms live) and **1.485e-10** (degenerate-eigenvalue edge rig); the apex tangent is rank 0 by construction so there is no FD gate on it; the iteration contrast is an OUTCOME — `Backward_Euler` with `Continuum` AND with `Secant` fails to converge (`analyze -> -3`) on the HB oedometric deck where `Closest_Point`/`Algorithmic` converges in **16**. gate 3 — nothing to gate: the family is registered perfectly plastic and the all-IVs-inert requirement is folded into the family marker at compile time. gate 4 — the 23-deck / 282-row `Backward_Euler` baseline still byte-identical; step refinement replaces P2's exact step independence, because a CURVED surface cannot be step independent — N = 1 IS the oracle (1.448e-15) and the sequence converges monotonically (**2.594e-3 → 1.413e-3 → 5.885e-4 → 1.260e-4** against N = 160) to a limit 0.26 % away. gate 6 — all seven deck-reachable MIXED Hoek-Brown pairings refused, each checked in BOTH directions; `strict_convergence` byte-inert (0.0); hydrostatic and near-degenerate tension states commit a finite vertex, not NaN. **Support 22 → 23 of 46:** `HoekBrown_YF` is registered in 7 specializations and `HoekBrown_PF` in 6, and in exactly ONE are both of the family; the other eleven stay refused, all pinned by `static_assert` in the g++ pre-flight. **The P3 decision, implemented and recorded:** `HoekBrown_PF::g` is evaluated in the un-negated frame and collapses to a **Tresca** potential (`HB_mb_psi` inert, exactly zero dilatancy, 32.86° off the normal at `mb_psi == mb`, no apex return past the tensile corner). `Closest_Point` uses a frame-consistent Hoek-Brown potential in its OWN path; the shipped `g` is untouched (D1) and the gap is PINNED in both directions — plastic `eps_vol` **+2.208111e-05** (the oracle's own +2.208e-05) against BE's **+2.09e-13**, a **5.910e-02** relative stress gap = **26.20 %** of the strength scale. Fixing `g` is the owner's separate PR; this WP is its warrant. **Also found:** `HB_sigma_ci` is not a parameter — the name is `HB_sigci` — so P1's and P2's Hoek-Brown refusal assertions were passing on a MISSING PARAMETER and never exercised the family gate; both corrected and inverted here. Full battery **204 passed, 2 skipped (pre-existing), 0 failed**. **gate 5 (mutation, [[reviews/adr97_p3_mutation]]):** dropping the `dl * D3 dm/dy` CURVATURE term from the Jacobian on a scratch build (one line in `hb_assemble`; `dm/dy` has exactly ONE non-zero entry and that entry IS the meridian's curvature) makes the free-node FACE rig stop converging entirely and takes the load-driven edge rig's FD error from **1.485e-10 to 8.309e-03** (5.6e+07x), while the local Newton goes 4/4/3/5 -> **7/8/7/7** iterations; the COMMITTED STRESS is unchanged (1.6e-15 ... 1.9e-13 against the oracle) because the residual is untouched, which is what makes the 21 survivors interpretable, and the P1/P2/P4/P6 + ADR-94 HB suites all still pass (86/86) because the mutation is confined to the Hoek-Brown tangent. 11 killed, 21 survivors. Reverted, rebuilt, re-measured identical. Full report: [[reviews/adr97_p3_report]].
## See also

[[94_asdplastic_review_plan]] · [[reviews/adr94_verdict]] (§1 M2/M3/M9, §6 D2, §7) ·
[[84_ladruno_mc_tension_cutoff_adr]] (the `special_return` hook and its SR tangent) ·
[[87_ladruno_depth_with_width_adr]] · [[LEDGER_vanilla_files]] ·
[[LEDGER_implementations]] · [[LEDGER_quirks]] (the ASDPlasticMaterial3D entries) ·
[[40_ladruno_performance_adr]] (P6's measurement harness)
