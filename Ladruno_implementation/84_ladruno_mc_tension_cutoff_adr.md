---
title: "ADR 84 — Mohr-Coulomb with tension cutoff for ASDPlasticMaterial3D"
project: Ladruno
type: ADR
status: P0 shipped; P2(a) shipped; P3 shipped (tangent-operator policy)
priority: high
owner: nmora
related: ["Cerro Lindo ADR-0005", "[[LEDGER_implementations]]", "[[LEDGER_vanilla_files]]", "[[Ladruno_materials_guide]]"]
tags: [adr, material, plasticity, geomaterial, asdplastic]
updated: 2026-08-13
---

# ADR 84 — MohrCoulombTensionCutoff for ASDPlasticMaterial3D

## 1. What / Why

The Cerro Lindo SSI program (consumer ADR `Cerro Lindo ADR-0005`, validation-ladder
rung M3) measured two defects in an ASDPlasticMaterial3D Mohr-Coulomb geomaterial,
per Gauss point, on a horseshoe-cavity + EDZ deconfinement model:

1. **Missing tensile cap.** The MC envelope extends into the tensile quadrant with
   intercept `c·cotφ` — ×11 the Hoek-Brown-derived tensile strength for their host
   rock (275 vs 24.7 kPa), unbounded relative to a fractured EDZ (~0). 498 EDZ
   exceedances at λ ≥ 0.95, max σ₁ = +282 kPa, clustered at the wall-foot corner
   where the M4+ rungs put the arch — an unconservative demand-transfer error.
2. **Inadmissible accepted states.** 20 GPs measurably OUTSIDE the yield surface
   (f/2c·cosφ = +0.0299): the Backward_Euler return map accepts f > 0 states near
   the MC apex.

This ADR adds a registered composite `MohrCoulombTensionCutoff_YF` +
`MohrCoulombTensionCutoff_PF` pairing with `TC_min_stress` as input, and the
integrator mechanism that returns robustly at the cutoff face/edge/corner/apex.

## 2. Root cause of defect 2 (verified in source, upstream code)

- `ASDPlasticMaterial3D.h` `Backward_Euler`: the apex block (`:2076` ff.) detects
  the apex region but its 60-line body is **entirely commented out** (references a
  nonexistent `pf.apex_flow_direction`); control falls through to the scalar
  Newton, where near the apex `n → sinφ/3·δ` is parallel to the volumetric part of
  `E·m` → `SINGULAR LOCAL TANGENT - FAILING!` or a stalled iteration.
- The Newton loop **accepts on exhaustion**: it falls out of `max_iter` and
  `return 0` with no convergence check — the persistence mechanism for f > 0.
- A second latent defect: the elastic early-exit fires whenever f merely
  *decreased* by more than `tol_yf` across the step. From admissible committed
  states it cannot accept f > 0, but from an already-inadmissible commit it can
  perpetuate one. Both upstream defects affect **all** ASDP materials and were
  **deferred out of P0** (changing them alters behavior for VonMises/DP users).
  **Both are now fixed opt-in under `strict_convergence` — see §6c (P2(a),
  shipped 2026-08-13)** and the two matching `LEDGER_quirks.md` entries.

## 3. Decision: composite YF + opt-in `special_return` hook, exact principal-space returns

### Alternatives considered

- **(A) Composite max-YF + revived apex hook.** The legacy `APEX_STRESS` signature
  takes no trial stress, so it cannot express predictor-dependent corner targets.
  Rejected as-is; its one-shot-projection *shape* survives in the new hook.
- **(B) Full exact principal-stress return (Clausen/Damkilde/Andersen) as an
  alternate integration branch.** Heaviest vanilla footprint; the MC-face return
  must stay on the existing scalar-Newton path anyway to keep existing-material
  behavior bit-identical. Rejected in pure form; its closed-form region returns
  survive as the targets inside the hook.
- **(C) Sequential/hierarchical returns.** No home: the framework is strictly
  single-surface (one yf object, scalar Φ, scalar λ in every scheme). Rejected in
  pure form; its escalation logic survives inside the hook.
- **Framework multi-surface machinery.** Far larger change to an upstream
  framework we do not own; not needed for a two-plane composite. Rejected.

**Chosen: hybrid.** One new YF class computing `f = max(f_MC, f_TC)` with
branch-consistent gradients, plus ONE new integrator mechanism: an opt-in trait
`yf_has_special_return` + `SPECIAL_RETURN` member that the Backward_Euler calls
after the elastic check and before the scalar Newton. The YF resolves every
cutoff-active trial state internally with exact principal-space returns.

### The geometric fact the design leans on

When `TC_min_stress < c·cotφ` (both Cerro Lindo cases), **the composite surface
has no MC apex**: the cutoff truncates it. Every trial state that used to reach
the degenerate apex neighbourhood has `f_TC > 0` and is intercepted by the hook,
landing on the Rankine plane whose face/edge/vertex returns are closed-form.
Conversely, any trial with `s3 ≤ T + tol` is provably safe for the untouched
scalar Newton: the MC return direction strictly decreases the largest principal
stress (isotropic E, coaxial return), so it cannot cross the cutoff.

### The escalation chain (zero-f>0 argument)

Classification happens ONCE from the trial stress (no per-iteration re-branching,
no corner chatter). For `f_TC > tol` trials, in order:

1. **Rankine face** (`s3' = T` exact — this is the uniaxial-tension acceptance);
2. **Rankine edge** (`s2' = s3' = T`, exact 2×2 in Lamé constants);
3. **MC∩TC corner** — Koiter two-multiplier return, directions frozen at the
   trial state (planes in principal space), residuals driven on the *framework's
   smoothed* f_MC and f_TC so the accepted point is admissible w.r.t. what the
   material actually evaluates. Non-associated MC flow (ψ) retained; associated
   Rankine flow on the cutoff.
   - `λ_TC < 0` ⇒ the true return is MC-only ⇒ **frozen-direction 1D MC Newton
     in-hook**;
3b. **Compound corner MC ∩ TC(s3) ∩ TC(s2)** — the triaxial-extension corner
   line of the composite, target CLOSED FORM: `(s1*, T, T)` with
   `s1* = (T(1+sinφ) − 2c·cosφ)/(1−sinφ)` (the smoothed f_MC collapses to the
   textbook form at θ=+30°, so the point is exact). This is **the attractor of
   symmetric lateral-tension drives** (s2 = s3): without it those states
   collapse to the vertex, losing the entire MC confinement (measured:
   −235.25 kPa of s1 replaced by +24.7). Cone membership is tested over the
   FULL four-direction fan {both MC faces meeting at the ridge, both Rankine
   flows} — testing a single MC face makes symmetric drives infeasible
   (μ-split condition `b·μΣ ≤ dep₂+dep₃`).
3c. **MC-dominant fall-through**: if every exact cutoff-feature return rejects
   the trial AND `f_MC ≥ f_TC`, the plastic increment is infeasible for every
   cutoff-feature cone because the true return lies on the MC surface BELOW
   the cutoff (TE-ridge returns from slightly-tensile lateral overshoots).
   Hand to the scalar Newton, whose MC-branch gradient treats it exactly as
   plain MC would; a converged Newton satisfies `max(f_MC, f_TC) ≤ tol` by
   construction. (Without this rung, transient global-iterate overshoots on a
   uniaxial-compression path latched the vertex and the global Newton
   limit-cycled at a fixed norm — measured, step-reproducible.)
4. **Apex/vertex** `T_eff·δ`, `T_eff = min(T, c·cotφ)` — always feasible by
   construction (`f_MC(T_eff·δ) ≤ 0` iff `T_eff ≤ c·cotφ`; `f_TC = T_eff − T ≤ 0`).
   Reached only by TC-dominant trials none of whose exact returns validated.

Every candidate is validated against both yield values to `tol_yf` before
acceptance; the chain terminates at the always-feasible apex. Worst case is a
small conservative stress drop (logged, first 5 events), never an accepted f > 0
state. **Structurally, zero inadmissible commits through this path.**

**Tangents:** all hook tangents are Koiter-consistent for their active set
(rank-1 face, rank-2 edge/corner, rank-3 compound/vertex).

> [!warning] P0 then **Secant-blended them inside the YF**, `(E + D)/2`, on the
> reasoning that raw rank-deficient tangents make the global stiffness
> near-singular. That was wrong in a way that shipped: it also **silently
> overrode the material's configured `tangent_type`** on every hook-resolved
> Gauss point, and at the corner it fabricated stiffness in the direction whose
> true tangent is ~0. It is the measured cause of the Cerro Lindo M5 stall.
> **P3 moves the policy to the integrator — see §9.**

### Cutoff above the apex (`T ≥ c·cotφ`)

The cutoff cannot truncate the MC cone: the material degenerates to plain MC.
`T_eff` clamps to the apex, a one-time warning is printed, and the hook takes
over exactly one region — the plain-MC apex cone — with a one-shot projection to
`c·cotφ·δ` (incidentally repairing the plain-MC apex defect for this material).

### Bit-identical strategy

- `f_MC` and the MC gradient branch are **verbatim copies** of `MohrCoulomb_YF`'s
  arithmetic; `f = max(f_MC, f_TC)` returns the winner's bits; the hook
  disengages (`return false`) whenever `s3 ≤ T + tol`. Confined-compression paths
  therefore reproduce plain-MC results bit-for-bit.
- Existing materials: every existing YF has `yf_has_special_return = false`, so
  the new integrator block compiles away — behavior is bit-identical by
  construction.

### Registration (explicitly NOT the cross product)

One explicit `COMBINED_MODELS` entry in the generator (mirroring
`STIFFSOIL_MODELS`): `LinearIsotropic3D_EL × MohrCoulombTensionCutoff_YF ×
MohrCoulombTensionCutoff_PF × BackStress<NullHardeningTensorFunction>`.
**+1 instantiation (45 → 46).** Pairings like `MohrCoulombTensionCutoff_YF ×
MohrCoulomb_PF` are deliberately not registered: dilatant MC shear flow on the
tension cap is mechanically wrong, and every extra pairing grows the
`OPS_AllASDPlasticMaterial3Ds.cpp` monster TU (see the C1060 quirk).

### Flow rule

ψ (non-associated) retained on the MC surface and at the corner; associated
(radial/Rankine) return on the cutoff — per the consumer's requirement. The PF's
cutoff branch matters only for tangent probes and future FE/RK support under the
default Backward_Euler; the hook computes its own multipliers.

### Deliberate non-replications / neighbours

- `MohrCoulomb_PF.h:83` scales cohesion by π/180 (benign there: additive constant
  in g, cancels in the differences). **Not replicated** in the new PF — there c
  enters the branch comparison and must be dimensionally right. Upstream file
  untouched.
- Standalone `TensionCutoff_YF` stays out of the generator (stub apex handlers,
  no PF, cross-product blowup). Deferred; the composite covers the requirement.
- The legacy dead apex block in `Backward_Euler` stays byte-identical (smaller
  upstream diff); the new hook is inserted after it.

## 4. Usage

```python
ops.nDMaterial("ASDPlasticMaterial3D", tag,
    "MohrCoulombTensionCutoff_YF", "MohrCoulombTensionCutoff_PF",
    "LinearIsotropic3D_EL", "BackStress(NullHardeningTensorFunction):",
    "Begin_Model_Parameters",
    "YoungsModulus", E, "PoissonsRatio", nu,
    "MC_phi", phi_deg, "MC_c", c, "MC_psi", psi_deg, "MC_ds", 0.0,
    "TC_min_stress", T, "MassDensity", rho,
    "End_Model_Parameters",
    "Begin_Internal_Variables", "BackStress", 0.,0.,0.,0.,0.,0.,
    "End_Internal_Variables")
```

> [!warning] ALL of `MC_phi, MC_c, MC_ds, MC_psi, TC_min_stress` are REQUIRED.
> ASDPlasticMaterial3D parameters are UNINITIALIZED unless set, and unknown names
> are a silent no-op. A forgotten `TC_min_stress` most likely reads 0.0 — a
> zero-tension rock — silently. (Candidate P2: opt-in required-parameter
> assertion in the factory.)

`updateParameter` ramps (the consumer's G-drive degradation clock) work with zero
extra wiring: the catch-all response ID 1000 routes any parameter name through
`setParameterByName`, and both YF and PF read parameters live per call. Tension
convention: TENSION-POSITIVE; `TC_min_stress` is the maximum allowed principal
stress (a tensile strength), sensible values `0 ≤ T < c·cotφ`.

Integration: documented for the default `Backward_Euler` (+ any tangent type).
Other schemes (FE/ME/RK45) run without the hook in P0 — see §6.

## 5. Acceptance / tests

Zone-A `tests/test_asdplastic_mctc.py` (t0m material-point driver, availability-
guard skip): uniaxial tension caps at exactly `TC_min_stress` (rel 1e-9, T=24.7
and T=0); confined-compression parity vs plain MC (bit-identical contract);
apex-regression on a hydrostatic-tension path (the test the M3 report says would
have caught the defect); corner probe across the MC∩TC intersection; unload-slope
tangent contract; updateParameter c/φ/T ramps; no-op cutoff (T > c·cotφ); and an
**f>0 audit** recomputing both yield values in numpy for every committed state of
every test — the M3 acceptance in miniature. Model-scale sign-off remains the
consumer's M3 re-run: zero f > 0 states, bounded EDZ σ₁, λ ≤ 0.9 unchanged.

Zone-A `tests/test_adr84_p3_confined_corner.py` (5 cases) adds the state the P0
battery never visited — **confined, non-associated MC shear on the cutoff face
and the MC∩TC corner at 4.86 MPa** — on a driver with **free DOFs**, without
which the material tangent is unobservable. See §9.

## 6. Phasing

- **P0 (this PR):** everything above under Backward_Euler.
- **P1 (measured):** wire the same hook into `Backward_Euler_LineSearch`
  identically if model-scale Newton statistics warrant; extend tangent-FD
  coverage to special-region states.
- **P2(a) — SHIPPED (this PR).** Upstream BE exhaustion-accept + f-decreasing
  elastic exit, fixed **opt-in** behind a new `strict_convergence` integration
  option (int, default 0). See §6c for the outcome.
- **P3 — SHIPPED (this PR).** Tangent-operator policy moved out of the YF into
  the integrator + `return_quality` / strict-mode refusal of the terminal vertex
  projection + the confined-corner Zone-A battery. Root-causes and closes the
  Cerro Lindo M5 stall. See §9.
- **P2 (still deferred):** (b) FE/ME/RK support via the same hook (their drift
  correction exists); (c) standalone TensionCutoff_YF (needs a real TC_PF +
  apex handlers); (d) apex fix for plain MohrCoulomb_YF; (e) required-parameter
  assertion.
- **P4 (new, demand-driven):** `Numerical_Algorithmic_*` differentiate
  `compute_local_stress()`, a simplified map that never calls `special_return`
  (nor the real `Backward_Euler`), so those operators are inconsistent with the
  actual response for **every** ASDP material. Framework-wide; out of P3's scope.

## 6b. Two defects found during P0 verification (both fixed/filed)

1. **Upstream NaN-surviving construction (FIXED here).** The ASDP constructor
   zeroed its Eigen members with `*= 0` — but they start as uninitialized heap
   storage, and garbage that decodes as NaN/Inf survives (`NaN*0 == NaN`).
   Fresh processes get OS-zeroed pages, so standalone probes always passed;
   pytest's churned heap recycles dirty blocks, so ~40% of suite runs got NaN
   in the shear slots of `CommitStress`. Every yf comparison on NaN silently
   goes false: plain MC trips its Newton NaN guard (analyze fails), while the
   MCTC escalation chain would land on the apex and return a *clean-looking*
   `T_eff·δ` on a compression path — precisely the garbage-into-plausible
   failure this feature exists to prevent. Fixes: constructor `setZero()`
   (all ASDP materials benefit), plus a finite-trial guard at the top of
   `special_return` (a non-finite trial falls through to the generic map's
   loud NaN check, never into the projection chain).
2. **`system("FullGeneral")` crashes the process on zero free equations**
   (every DOF fixed or sp-prescribed) — material-independent (measured with
   ElasticIsotropic, plain MC, MCTC), `UmfPack` returns rc=0 on the identical
   model. Pre-existing core defect, out of scope here: filed as its own task +
   quirks entry; the Zone-A driver uses UmfPack.

## 6c. P2(a) outcome — `strict_convergence` (shipped 2026-08-13)

Both §2 defects are fixed behind one opt-in integration option,
`strict_convergence` (int, default 0, parsed in `Begin_Integration_Options`,
stored per material tag in `INT_OPT_strict_convergence`). With it ON, (a) the
scalar-Newton loop tracks a `be_converged` flag and, on exhaustion with
`|Phi| >= tol_yf`, prints an `opserr` line (tag, final `|Phi|`, tol) and returns
-1 instead of falling through to `ComputeTangentStiffness(); return 0;`; and (b)
the elastic early-exit's f-decreasing disjunct additionally requires
`yf_val_end <= tol_yf`, so an inadmissible commit is corrected by the plastic
corrector rather than carried forward (elastic unloading, `f_end < 0`, still
exits). Default 0 is byte-identical to upstream, deliberately: turning it on
changes convergence behaviour for every existing VonMises/DP/MC deck. The P0
`special_return` hook sits upstream of the Newton loop and is unaffected in
both modes. **Measured:** the whole shipped ASDP battery
(`test_asdplastic_mctc.py` + `test_asdplastic_response_tags.py`) passes
unchanged with the flag off (3x, for the churned-heap NaN trap) **and** passes
unchanged with `strict_convergence 1` force-injected into all 19 of its
ASDPlasticMaterial3D definitions — zero behaviour differences on healthy paths.
New Zone-A battery `tests/test_adr84_p2a_strict_convergence.py` (6 cases).

Three findings came out of the measurement, all pinned by that battery:

1. **The apex path is not the reproducer.** Plain `MohrCoulomb_YF` pulled 1.5x
   past its hydrostatic apex — the obvious candidate, and the ADR's own
   documentation run — converges cleanly on this build: 30/30 steps, worst
   committed `f_MC = 2.8e-14`, flag-on bit-identical. What reproduces the
   exhaustion-accept deterministically is starving `n_max_iterations` on a
   plainly plastic leg, which is the defect stated precisely rather than a
   contrivance. Measured at `n_max_iterations 2`: flag off completes 20/20 steps
   *all reporting success* with **all 20 committed states inadmissible** and
   worst `f_MC = 77.6 kPa` against a 2.8e-1 tolerance (276x); flag on refuses at
   step 1 and commits nothing.
2. **`stdBrick` swallows the refusal — a THIRD silent accept, in a vanilla
   element.** `Brick::update()` writes `success = ...->setTrialStrain(strain);`
   then unconditionally `return 0;` — assigned, never read. So the material
   refuses, prints, returns -1, and the analysis still reports success.
   `TenNodeTetrahedron` accumulates and returns the sum (already fixed in the
   fork — TIMs report item 8), so the contract tests use the tet host. NOT fixed
   here on purpose: `return success` is an unconditional behaviour change for
   every stdBrick + every material, exactly the blast radius the opt-in flag
   exists to avoid. Filed as its own follow-up; pinned by
   `test_stdbrick_swallows_the_refusal` so fixing it surfaces as a loud,
   informative failure.
3. **Strict mode gates every global-Newton TRIAL, not only the commit — so it
   needs a stress-scaled `f_absolute_tol`.** A step's first global iterate can be
   far from the solution and the local return map on that large trial increment
   may legitimately exhaust. Measured on the tet rig (max |sigma| ~ 2.8e5 kPa,
   `n_max_iterations 100`): at `f_absolute_tol` 1e-6 (the default), 1e-4 and
   1e-2 strict mode refuses step 1 of a run that is healthy flag-off; at 1e-1 —
   a ~3.6e-7 *relative* demand — flag-on is bit-identical to flag-off. **Operating
   guidance of record: scale `f_absolute_tol` to your model's stress magnitude
   before turning `strict_convergence` on.** The absolute-tolerance design is
   upstream's and is left as-is here; making the yield tolerance relative is a
   separate decision with its own blast radius.

## 7. Risks

- **Corner 2×2 with frozen directions vs DP-smoothed curvature:** bounded
  iterations, escalation to the always-feasible apex; worst case conservative,
  logged. No f > 0 escape.
- **Eigenvalue coalescence at the cutoff:** face return with s2 ≈ s3 escalates to
  the edge return, symmetric in the degenerate subspace (unique result); FD
  fallback covers gradient probes.
- **ψ = 0 (EDZ) / extreme ν:** corner Jacobian dominated by
  `n_TC·E·m_TC = λ_L + 2G > 0` — nonsingular.
- **Monster TU:** +1 instantiation is trivial; the C1060 serial-first recipe in
  `LEDGER_quirks.md` still applies.

## 8. Ledger / obligations

- `LEDGER_implementations.md`: one feature row (new YF/PF headers + test + this
  ADR). No new class tag — the ASDP template family shares `ND_TAG_
  ASDPlasticMaterial3D = 10000`.
- `LEDGER_vanilla_files.md`: rows for `YieldFunctionBase.h`,
  `ASDPlasticMaterial3D.h`, `AllYieldFunctions.h`, `AllPlasticFlowDirections.h`,
  `gen_ASD_material_definitions_CPP.py`, `ASD_material_definitions.cpp` — all
  marked `// Ladruno (ADR-84 P0)`.
- `LEDGER_quirks.md`: entries for the BE exhaustion-accept and f-decreasing
  elastic-exit upstream defects, and for the silent-uninitialized-parameter trap.
- Banner line in `Ladruno_scripts/banner_features.txt` + `patch_banner.py`.
- Upstream campaign: candidate package alongside 1.5 (Hoek-Brown/StiffSoil) —
  the composite YF/PF are clean additive headers on the stock ASDP framework;
  the hook block is the only integrator edit and is `if constexpr`-gated.

## 9. P3 — the confined-shear stall (Cerro Lindo M5), root cause and fix

**Report (13 Aug 2026).** The M3/M4 model (2-D plane strain, `LadrunoLST`/
`LadrunoQuad`, host rock c = 1014 kPa / φ = 45.95° / ψ = 11.49°, EDZ c = 100 /
φ = 20° / ψ = 5°, in-situ σ₀ = 6160 kPa) runs to λ = 1 with plain MC. With the
composite substituted and nothing else changed it stalls in **stage S1 — pure
continuum** — at λ ≈ 0.45-0.50: 3 819 PARDISO pivot-perturbation warnings, and
Newton sitting at a ~9e-5 m displacement increment against a ~500 kN residual.
Three settings, all stall: shipped defaults (λ = 0.500); `strict_convergence` +
`f_absolute_tol` 1e-1 + 200 iterations (λ = 0.450); `tangent_type = Continuum`
(λ = 0.475). Their reading: the corner between the shear surface and the cutoff
is active, and the consistent tangent there is near-singular.

**The reading was right about the region and wrong about the mechanism, and the
"all three settings behave the same" clue was the actual tell.**

### 9.1 What was measured

A numpy replica of the shipped algorithm (`special_return` chain + the scalar
Newton, transcribed from the merged C++) was used to reach the states a global
Newton actually presents, then compared against numerical differentiation of
the material's own return map.

1. **Corner chatter is NOT the mechanism.** Over a 45×45×3 sweep of trial states
   at 1-6 MPa confinement, the scalar Newton failed to converge on 17 trials —
   **all with zero branch flips** (every iterate on the MC branch). The
   `f_TC ≥ f_MC` branch switch does not oscillate; the escalation chain's
   once-only classification holds.
2. **The corner is a rigid attractor.** Reached from an isotropic in-situ state,
   the return pins the stress at `(-4862.24, -4289.39, +24.70)` **to the last
   digit for a 20× range of strain increment**. The true algorithmic tangent
   there has `min|eig| = 9.4e-4`, `cond = 6.8e9`.
3. **The shipped tangent fabricates stiffness in exactly that direction.**

   | operator | ‖D‖ | rel. err vs truth | cond | min\|eig\| |
   |---|---|---|---|---|
   | numerical (truth) | 5.73e6 | — | 6.8e9 | **9.4e-4** |
   | raw Koiter rank-2 | 6.08e6 | **0.121** | 9.0e14 | 5.8e-9 |
   | **secant blend (P0, shipped)** | 8.14e6 | **0.868** | 3.2 | **2.0e6** |
   | elastic | 1.20e7 | 1.724 | 5.0 | 2.0e6 |

   The blend is 87% wrong and, worse, **well-conditioned where the truth is
   rank-deficient**. Plain MC on the same path hands out `cond ~1e16` and
   converges. A global Newton told "push here and stress rises by 2e6·Δε",
   whose stress then does not move at all, produces precisely the reported
   signature: a vanishing increment against a residual that never falls.
4. **`tangent_type` was inert.** `special_return` wrote `stiffness_return`
   already blended and the integrator assigned `Stiffness = stiff_sr`
   unconditionally, so **every hook-resolved Gauss point got the blend no matter
   what the deck asked for**. That is why the consumer's three settings all
   behaved alike — the knob they reached for could not reach the Gauss points
   that needed it.
5. **Not corner-specific.** The Rankine FACE rung is blended too; on the P3
   regression path the face steps carry the same penalty. The M5 model has far
   more face points than corner points, which fits a stall in *pure continuum*
   at λ ≈ 0.45 rather than at a small corner cluster.

### 9.2 Fix

The YF no longer applies a tangent policy — `stiffness_return` is the **raw
active-set (Koiter) tangent** — and `Backward_Euler` applies the material's
configured operator to it, exactly as the generic path does:
`Elastic → E`; `Continuum`/`Algorithmic` → raw; `Secant` (default) → `(E+D)/2`.
`Numerical_Algorithmic_*` deliberately fall back to the raw active-set tangent:
they differentiate `compute_local_stress()`, a simplified map that never calls
this hook, so they would be strictly worse here (a pre-existing framework
inconsistency, not addressed here).

**Default is `Secant`, so this is byte-identical unless the deck asks
otherwise.** The whole shipped ADR-84 battery passes unchanged (20/20).

`special_return` also gained an out-param `return_quality`
(`SR_QUALITY_EXACT` / `SR_QUALITY_FALLBACK`). The terminal Stage-4 vertex
projection — the one inexact rung, which replaces a confined deviatoric state
with `T_eff·δ` — now reports `FALLBACK`, and under `strict_convergence` the
integrator **refuses the step loudly** instead of committing it silently. That
is the M5 report's request #3.

### 9.3 Measured effect (`tests/test_adr84_p3_confined_corner.py`)

Confined driver, host-rock parameters, corner at **4.86 MPa confinement**
(`s1 = -4862.2`, `s2 = 0`, `s3 = T = 24.7`) — inside the σ₃ ≈ 1-6 MPa band the
report asked for; 45/45 steps with the cutoff active, 15 on the corner:

| `tangent_type` | total Newton iterations | per working step |
|---|---|---|
| Elastic | 294 | 9 |
| **Secant (default)** | 232 | 7 |
| **Continuum** | **76** | **2** |

**3.5× fewer iterations** on the 31 working steps (62 vs 218).

### 9.4 Operating guidance for the consumer

Set **`tangent_type Continuum`** on the MCTC materials. Before P3 that setting
was inert on the hook-resolved Gauss points; it is now the raw active-set
tangent, which is the correct Newton operator on the cutoff face and at the
corner. Keep `strict_convergence` on with a stress-scaled `f_absolute_tol` (§6c
finding 3) so a vertex-fallback or a non-converged return map refuses loudly
rather than committing a quietly wrong state.

### 9.5 Why P0 could ship this

Two gaps in the P0 battery, both closed by the P3 module:

- it exercised the composite in **uniaxial tension** only — never confined,
  non-associated MC shear at the corner;
- its material-point driver **prescribes every DOF**, so there are zero free
  equations, the global Newton never solves anything, and the tangent the
  material hands the assembler is *unobservable*. `test_tangent_fd` checks the
  tangent by finite differences of the material response, which the blend
  passes — it is a self-consistent operator, just not the right one.

The P3 driver leaves the z-faces free precisely so the tangent shows up in
`testIter()`. **Any future tangent work on this material must use a driver with
free DOFs.**

## See also

- Clausen, J., Damkilde, L., Andersen, L. (2006/2007) — exact MC returns in
  principal stress space, incl. tension cutoff (return-region taxonomy).
- Koiter, W.T. (1953) — multi-surface plastic flow at corners.
- `Models/SSI FEM/NOTE_M3_RevA.md` G3 (consumer measurement),
  `NOTE_M0-M2_RevA.md` trap 9 (the declared-but-inactive cutoff discovery).
