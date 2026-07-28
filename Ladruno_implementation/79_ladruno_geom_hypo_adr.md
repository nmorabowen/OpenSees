# ADR 79 — `-geom hypo`: hypoelastic rate-form updated-Lagrangian geometry method

**Status:** P0+P1 shipped (dry solid); P2 (u–p) staged, blocked on ADR-78 (PR #677) merge
**Element:** `LadrunoBrick` (ELE 33002) first; `LadrunoUP` (ELE 33017) in P2
**Seam extended:** `SolidTransformation` gains `METHOD_HYPO = 3` (marker only)
**Kernel:** `SRC/element/solidTransformation/LadrunoHypoKernel.h` (header-only, OpenSees-free)

---

## 1. Why

ADR 78 opened `-geom corot` on the u–p family: large ROTATION + geometric
stiffness, small deformational strain. The remaining gap is genuine LARGE
STRAIN for soils. `-geom finite` (the `setTrialF(F)` path) is blocked on the
material side: the `LogStrainNDMaterial` lifting adaptor states its own limit —
*"v1 assumes a LINEAR-elastic inner law (so Cᵉ = inv(initial tangent) is
constant and the εᵉ recovery is exact)"* — and the target soil model
`PressureDependMultiYield` is pressure-dependent hypoelastic (G ∝ p′^0.5),
which violates that assumption directly. No finite-strain soil model exists in
the fork or upstream, and writing one is a multi-year constitutive program.

The industry answer — what Abaqus NLGEOM does for every rate-form
(hypoelastic) material (Abaqus Theory Guide §1.4.3 strain-increment
integration, §1.5.3 objective rates, §1.5.4 state storage; via the
`abaqus-theory-foundations` skill) — is an **objective-rate updated-Lagrangian
integration**: compute an objective (rigid-rotation-free) strain INCREMENT per
step from midpoint-configuration gradients (Hughes & Winget 1980), hand it to
the SAME small-strain constitutive routine, keep the material's tensorial
state consistent across rotation, and assemble on the CURRENT configuration
with the initial-stress geometric term. PDMY-class soil laws are already
rate-form — they consume strain increments and produce stress increments —
so this route gives large strain WITH the existing soil library, unchanged.

---

## 2. The design decision: wrapper vs seam vs element lane

Three candidate homes for the hypo path were on the table. The deciding
constraint fell out of a fact about OpenSees materials, not about elements:

**The classic Hughes–Winget prescription is unimplementable against a stock
`NDMaterial`.** HW (and Abaqus's UMAT contract) rotate the stored stress AND
every stored tensorial internal variable (backstress, surface centers) forward
by the incremental rotation ΔR at the top of each step. That works when the
host code owns the state arrays (Abaqus stores STRESS/STATEV and hands them to
the routine). An OpenSees `NDMaterial` owns its state privately —
`PressureDependMultiYield`'s nested-surface centers, a J2's backstress — and
exposes no "rotate your tensors by ΔR" entry point. Adding one per material
would mean touching the whole library, which is exactly what this ADR exists
to avoid.

**Resolution: use the unrotated-configuration (Green–Naghdi) variant of the
same objective integration** (Flanagan–Taylor lineage; the PRONTO/Sierra
approach; Abaqus itself uses the Green–Naghdi rate for all Explicit VUMAT and
structural elements — Theory Guide Table 1.5.3-1). Instead of co-rotating the
material state forward each step, DE-ROTATE the objective strain increment
into a fixed material frame:

    Δε_mat = Rᵀ_mid · Δε_spatial · R_mid ,     R = polar(F)

and let the material live its entire life in that fixed unrotated frame —
its stress, its backstress, its surface centers never need rotating, because
its input never rotates. The stress is pushed forward `σ = R σ_mat Rᵀ` only at
assembly time. Objectivity class is identical (Jaumann and Green–Naghdi differ
only under combined finite rotation + finite shear; GN's simple-shear response
is monotone where Jaumann's famously oscillates — a strictly better ceiling).
And the fixed-material-frame convention is exactly the one `-geom corot`
already established (`SolidTransformationCorot.h` header notes: the material
and its backstress live in a fixed reference frame), so the kinematic-hardening
objectivity test pattern carries over unchanged.

With that settled, the three candidates:

1. **A fourth `SolidTransformation` method (seam-resident logic).** Rejected.
   The seam's LOCKED division of labour gives the transformation exactly four
   jobs (strain measure, localizeDisp, globalizeForce/Stiff, recorder frame) on
   whole-element vectors. The hypo path needs PER-GAUSS-POINT machinery —
   midpoint gradients, per-GP polar rotation, per-GP accumulated strain — that
   the seam interface cannot see (its `update()` receives only nodal
   coordinates). Cramming per-GP state into the wrapper would break the
   "formulation owns B" contract and make the wrapper element-shape-aware.
   Precedent: `-geom finite` hit the same wall and became an element lane with
   `SolidTransformationFinite` as a pure MARKER (all seam methods identity).

2. **A per-GP material wrapper (`HypoIncrementalNDMaterial`).** Rejected. The
   wrapper would need Δε and R per step — data outside the `NDMaterial`
   interface — so the element must `dynamic_cast` and call custom methods
   anyway; the wrapper buys no decoupling. It costs a user-facing extra object
   (`nDMaterial HypoWrap 2 1` around every soil material — the LogStrain
   friction, again), a full recorder/response passthrough layer (the SEAM-1
   lesson), and an ownership split of the per-GP state. With the unrotated-
   frame resolution above, there is nothing material-side left to wrap: the
   inner material is driven by plain `setTrialStrain()` with an accumulated
   feed — the element can do that directly.

3. **An element-internal lane (the `-geom finite` pattern).** **Chosen.**
   `SolidTransformationHypo` is a marker (METHOD_HYPO = 3, all seam methods
   identity, `getStrainMeasure() == SmallStrain` — honest: the material IS
   driven via `setTrialStrain`), parse/serialize rides the existing method-id
   plumbing, and the element implements `updateHypo()` +
   `formResidAndTangentHypo()` beside `updateFinite()` +
   `formResidAndTangentFinite()`, reusing the UL assembly skeleton (spatial
   gradients via F⁻¹, current-config dv = J·dV, initial-stress geometric term).
   The per-GP math lives in a header-only, OpenSees-free kernel
   (`LadrunoHypoKernel.h`, the `LadrunoUPKernel.h` discipline) so P2 reuses it
   from `LadrunoUP` without touching the brick.

**User experience result (the point of the whole exercise):**

    nDMaterial PressureDependMultiYield 1 ...
    element LadrunoBrick 1  n1..n8  1 -geom hypo

— no adaptor object, no material change, large strain.

---

## 3. The algorithm (per Gauss point, per evaluation)

State: committed accumulated feed strain `ε_feed_n` (6 Voigt components/GP —
the ONLY persistent hypo state) plus whatever the material commits internally.
Committed nodal displacement `u_n` comes from `Node::getDisp()` (committed);
trial `u_{n+1}` from `getTrialDisp()`. Every Newton iteration recomputes from
committed state — no accumulation-within-iteration hazard.

1. **Deformation gradients** from reference shape gradients `∂N/∂X`:
   `F_n = I + Σ u_n,a ⊗ ∂N_a/∂X`, `F_{n+1}` likewise, `F_mid = ½(F_n + F_{n+1})`
   (Hughes–Winget's midpoint configuration `x̄ = ½(x_n + x_{n+1})`).
2. **Midpoint spatial gradients** `∂N_a/∂x̄ = (∂N_a/∂X)·F_mid⁻¹`; incremental
   displacement gradient `G = Σ Δu_a ⊗ ∂N_a/∂x̄`, `Δu = u_{n+1} − u_n`.
3. **Objective strain increment** `Δε_spatial = sym(G)`. HW's midpoint makes
   this EXACTLY zero for an arbitrarily large rigid-rotation increment (the
   sym part of `2·tan(Δθ/2)·ŵ` vanishes identically — verified in the kernel
   oracle), not merely O(Δθ²); the O(Δθ³/step) error quoted for HW applies to
   combined stretch+rotation paths.
4. **Midpoint rotation** `R_mid = polar(F_mid)` (Higham iteration, the corot
   wrapper's algorithm). De-rotate: `Δε_mat = R_midᵀ · Δε_spatial · R_mid`.
   Using the TOTAL polar rotation (not a ΔR product chain) keeps R stateless,
   drift-free, and restart-safe.
5. **Accumulate and feed:** `ε_feed = ε_feed_n + Δε_mat` (Voigt, engineering
   shear), then `setTrialStrain(ε_feed)`. The material perceives exactly the
   increment `Δε_mat` relative to its committed state — the Green–Naghdi-rate
   constitutive integration, by construction, on an UNCHANGED material.
6. **Assembly** (current configuration, `R_{n+1} = polar(F_{n+1})`):
   - push forward `σ = R σ_mat Rᵀ`;
   - `f_a = ∫ σ·∇ₓN_a dv`, `dv = J·dV`, spatial gradients via `F_{n+1}⁻¹`;
   - `K = ∫ BᵀcB dv + ∫ GᵀΣG dv` with `c` = the material 6×6 tangent
     bond-rotated by `R_{n+1}` and the standard symmetric initial-stress term
     `(K_geo)_{(a,i)(b,k)} = δ_ik ∫ ∂Nₐ/∂x_j σ_jl ∂N_b/∂x_l dv`;
   - body force stays a dead load per REFERENCE volume (the finite-lane
     convention, consistent with reference-config mass).

`commitState()` copies the trial feed to `ε_feed_n`; revert restores it;
`sendSelf` ships it (§8).

### 3.1 Conventions, stated honestly

- **Rate law on Cauchy stress** (`σ̇_GN = C:d`), the FLAC/geotech convention.
  No `J`/work-conjugacy (Kirchhoff) correction: at 20 % volumetric strain the
  difference from a Kirchhoff-based law is O(10 %) — that is the accepted
  hypoelastic trade, and it is what PDMY's calibration assumes anyway (soil
  moduli are fitted to true-stress/engineering-ish data at exactly these
  strains). The P1 gate vs `-geom finite`+LogStrain MEASURES this 1/J factor
  and asserts it, turning the discrepancy into a verified property.
- **Hypoelastic elasticity is not hyperelastic**: closed strain cycles through
  large strain dissipate/generate O(ε³) energy. Documented ceiling; irrelevant
  for the driving application (soils are dissipative from the first increment).
- **Simple shear at very large γ**: the Green–Naghdi rate is monotone (no
  Jaumann oscillation) but stiffens vs the analytic hyperelastic reference for
  γ ≳ 1. Gate 5 pins the element against the numpy oracle integrating the SAME
  algorithm — element-vs-oracle is exact; oracle-vs-physics is the documented
  hypoelastic ceiling.
- **Tangent policy**: `∫BᵀcB + ∫GᵀΣG` omits the rate-convection terms of the
  exact GN consistent tangent (O(σ) vs O(C), unsymmetric). Same "residual
  exact, tangent dominant-term" policy as corot v2.0 / ADR 78 §3.3; K stays
  SYMMETRIC (a real solver advantage over F-bar finite). Newton converges to
  the correct equilibrium; only the asymptotic rate is affected.
- **Recorder frame**: `stress`/`strain` material responses report in the
  UNROTATED material frame (σ_mat, ε_feed) — the corot convention, and what
  the objectivity gates want to see invariant. `frameTimeVarying()` stays
  false in v1 (per-GP frames don't fit the element-level recorder contract);
  a global-frame push-forward channel is a demand-driven follow-up.

---

## 4. Scope and guards (v1)

- `-formulation std` only. bbar/F-bar-style mean dilatation in rate form is a
  known extension (apply the mean-dilatation split to `Δε` and the assembly B)
  but is NOT needed by the driving application: drained PDMY is not
  incompressibility-locked, and undrained saturation arrives via the u–p
  element in P2, where incompressibility is carried by the pressure FIELD, not
  by the skeleton operator. uri/ssp/eas: reserved (same status as corot).
- Requires a SMALL-STRAIN material — a `FiniteStrainNDMaterial` is rejected at
  parse time (the existing GEOM-1 converse guard already does this for any
  non-finite method).
- `-damp` reserved (mirror of the finite guard); `-bulkViscosity` stripped
  (existing non-linear-geom guard).
- 3D H8 lane only in v1 (the kernel is dimension-generic; 2D and Tet lanes are
  demand-driven).
- `-geom linear` stays bit-identical: every hypo branch is behind
  `isHypo()`, which is false unless METHOD_HYPO was parsed.

---

## 5. Staging

- **P0 — kernel + oracle (this PR):** `LadrunoHypoKernel.h` (midpoint
  gradients, Higham polar, Δε build/de-rotation, Voigt bond rotation of σ and
  C) + numpy oracle `tests/hypo_reference.py` + g++ cross-check
  `tests/test_hypo_kernel_cpp.py` (the `test_ladruno_up_kernel_cpp.py`
  pattern: oracle emits cases, C++ harness diffs to ≤1e-9, plus standalone
  invariants — rigid-increment Δε ≡ 0, polar orthogonality, rotation-of-
  identity-tangent invariance).
- **P1 — dry solid (this PR):** the hypo lane on `LadrunoBrick` (H8, std).
  Gates (tests/test_ladrunoBrick_hypo.py):
  1. rigid-rotation zero stress through large multi-step rotation (tight —
     §3 step 3 makes it near-identity, asserted at solver tolerance);
  2. `-geom linear` unchanged + `-geom hypo` ≡ linear at infinitesimal strain;
  3. agreement with `-geom corot` at moderate strain under rotation;
  4. uniaxial-strain elastic compression to 25 %: exact vs the analytic
     hypoelastic target `σ = (λ+2μ)·ln λ`, and the measured gap to
     `-geom finite`+LogStrain asserted ≈ the 1/J factor (§3.1);
  5. elastic simple shear to γ = 2 vs the numpy oracle (exact) — documents
     the GN signature (monotone, no Jaumann oscillation);
  6. J2 kinematic-hardening cyclic objectivity in the hypo path's OWN
     convention (fixed material frame — the corot test pattern).
- **P2 — the u–p element (follow-up PR, after #677 merges):** extend to
  `LadrunoUP` reusing the kernel: the ADR-78 machinery carries over (GN22/GN11
  incremental p-row — the chord-defect findings), the incremental volumetric
  term becomes `Δ ln J` (continuity on the true volume change), and the
  ADR-78 §3.7 deferrals finally open: porosity `n = 1 − (1−n₀)/J` and
  Kozeny–Carman `k(J)`. Gates: ADR-78 battery rerun, drained equivalence vs
  the P1 dry hypo brick (roller sides; `-b` conventions differ — UP is an
  acceleration ×ρ_sat, brick is force/volume; summed-reaction weight check),
  large-strain 1D consolidation vs a Gibson-family reference.
- **P3 — PDMY01 smoke** on the TIMs footing regime (`References/SFIM_model/`,
  bearing limit point) — the driving application, last.

---

## 6. Serialization

Per-GP persistent hypo state = committed `ε_feed_n` only (8 GP × 6 = 48
doubles). Shipped as an EXTRA guarded `Vector(48)` after the standard `dData`,
gated on the geometry method id in `idData(28)` (the EAS-alpha pattern —
non-hypo streams are byte-identical; eas+hypo cannot co-occur, parser-refused,
so the guarded-send order is deterministic). `R` is stateless (recomputed from
nodal positions), so nothing else ships. The idData(28) packing
(`+ 1000·methodID`) already accommodates id 3 without widening.

---

## 7. Files

| file | change |
|---|---|
| `SRC/element/solidTransformation/SolidTransformation.h` | `METHOD_HYPO = 3` |
| `SRC/element/solidTransformation/SolidTransformationHypo.h/.cpp` | marker method (identity seams) |
| `SRC/element/solidTransformation/SolidTransformationLinear.cpp` | `create()` case 3 |
| `SRC/element/solidTransformation/LadrunoHypoKernel.h` | P0 kernel (header-only) |
| `SRC/element/ladrunoBrick/LadrunoBrick.h/.cpp` | hypo lane: state, updateHypo, formResidAndTangentHypo, commit/revert/send/recv |
| `SRC/element/ladrunoBrick/OPS_LadrunoBrick.cpp` | `-geom hypo` parse + guards |
| `tests/hypo_reference.py` | numpy oracle (emits `tests/hypo_cases.txt`) |
| `tests/hypo_kernel_check.cpp` + `tests/test_hypo_kernel_cpp.py` | P0 g++ cross-check |
| `tests/test_ladrunoBrick_hypo.py` | P1 gates |

No upstream file is touched; no new class tag (the seam is tag-free by
design).

---

## 8. References

- Hughes, T.J.R. & Winget, J. (1980), "Finite rotation effects in numerical
  integration of rate constitutive equations arising in large-deformation
  analysis", IJNME 15.
- Flanagan, D.P. & Taylor, L.M. (1987), "An accurate numerical algorithm for
  stress integration with finite rotations", CMAME 62 (the unrotated-
  configuration formulation).
- Abaqus Theory Guide §1.4.3 (strain-increment integration), §1.5.3 (Jaumann
  vs Green–Naghdi, Table 1.5.3-1), §1.5.4 (state storage / ΔR co-rotation) —
  via the `abaqus-theory-foundations` skill.
- de Souza Neto, Perić & Owen (2008) §14 — the contrast case: the
  multiplicative/log-strain path (`-geom finite`), and why it needs a
  linear-elastic inner law in the fork's adaptor.
- ADR 09/10 (SolidTransformation seam, corot), ADR 78 (u–p corot, the
  incremental p-row findings that P2 inherits).
