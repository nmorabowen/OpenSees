---
title: Solid transformation wrapper — the element geometry-method layer
project: Ladruno
status: draft
priority: high
owner: nmora
tags:
  - implementation
  - element
  - solid
  - kinematics
  - corotational
  - finite-strain
  - contract
---

# Solid transformation wrapper — the element geometry-method layer

A reusable **geometry-method layer** for solid (continuum) elements that
decouples *kinematics* from the *constitutive + integration core*. It is the
solid-element analogue of OpenSees' beam-column `CrdTransf` / `CorotCrdTransf` —
which OpenSees has for **frames only**, not for solids. First consumer:
[[09_ladruno_brick]]. Not called "corotational" on purpose: corotational is just
**one** method this layer hosts; total-Lagrangian finite strain is another, and
linear (identity) is the third.

Companion docs: the element that consumes it [[09_ladruno_brick]]; the recorder
frame contract [[ladruno_element_contract]] §B (per-step `localAxes`); the
**material-side log-strain adaptor** (owned by the parallel material-wrapper
work — referenced here as *seam 3*, see *The boundary with the material wrapper*);
and the **staged-activation hook** [[staged_deformation_gradiend]] — a capture of a
per-GP birth deformation gradient `F₀` on top of the `finite` ledger here, so an
element appended mid-stage is born stress-free at the deformed geometry (`F_rel =
F·F₀⁻¹`). Orthogonal to the linear/corot/finite *method* axis; see also the SP/MP/
element reference-frame trichotomy in [[constraints_reference_position]].

> **The one idea.** An element that knows its reference nodal positions **X** and
> its nodal field **u** already has everything needed to compute *any* strain
> measure. The layer is therefore **not** new information — it is a discipline for
> *where* that computation lives, so the constitutive core stays oblivious to
> which kinematics wrap it. Linear ⊂ corotational ⊂ finite is a *method* axis,
> orthogonal to the element's *formulation* axis (`std`/`bbar`/`uri`/`eas`).

## What

A small interface, `SolidTransformation`, that sits between global DOFs and an
oblivious small-strain element kernel. It performs three jobs each
residual/tangent evaluation:

1. **Kinematics ledger (intrinsic):** from (**X**, **u**, ∇N) at each Gauss
   point, produce the strain quantity the core/material consumes. This step
   *exists in every method* — it is the element's core duty, not an add-on.
2. **De-rotate / localize (method-specific):** transform global nodal kinematics
   into the form the small-strain core expects.
3. **Re-rotate / globalize (method-specific):** map the core's internal force and
   stiffness back to global DOFs, **adding the geometric stiffness** the method
   implies.

Three concrete methods:

| Method | Handles | Strain to material | Adds K_geo | Material library |
|---|---|---|---|---|
| **`linear`** | small disp + small rot | engineering ε = B·u | no | full (small-strain) ✓ |
| **`corot`** | **large disp + large rot**, small strain | small ε in de-rotated frame | yes (from frame ∂R/∂u) | full (small-strain) ✓ |
| **`finite`** (spatial multiplicative) | large **strain** | F → εᵉ = ½ln **b**ᵉ (spatial/Eulerian log strain) | yes (spatial geometric) | small-strain library *promoted via* seam-3 log-strain adaptor (Kirchhoff **τ** → Cauchy **σ**) |

## Why

- **Large rotation without a rewrite.** Most "geometric nonlinearity" in
  earthquake/structural work is large-rotation / small-material-strain (P-Δ,
  member buckling, soft-story sway, pounding rotation). `corot` covers it while
  **reusing the entire small-strain material library** (feed the material small
  corotated strain, get Cauchy stress back) and **every formulation in the core**
  (`bbar`/`uri`/`eas` run unchanged inside the de-rotated frame).
- **Large strain when it's actually needed** (geotech large-deformation &
  post-failure, base-isolation elastomers, localization) — `finite` + the
  log-strain material adaptor promotes the existing small-strain soil/metal
  models to large strain (see *seam 3*).
- **No solid `CrdTransf` exists upstream.** OpenSees' transformation abstraction
  is beam-column-only. Solids that want geometric nonlinearity today must be
  rewritten per-element. This layer gives solids the same clean separation frames
  already enjoy.
- **Verifiability.** Each method is testable in isolation: `linear` ≡ the
  small-strain element; `corot` must give **zero stress under arbitrary rigid-body
  rotation** (the defining patch test); `finite` reduces to `corot` then `linear`
  as deformation → 0.

## Where

- **New code:**
  - `SRC/element/solidTransformation/SolidTransformation.h` — abstract interface.
  - `SolidTransformationLinear.{cpp,h}` — identity (the default; pure passthrough).
  - `SolidTransformationCorot.{cpp,h}` — polar-decomposition corotational + K_geo.
  - `SolidTransformationFinite.{cpp,h}` — total-Lagrangian (F, E, initial-stress
    K_geo); pairs with the seam-3 material adaptor.
- **Consumed by:** [[09_ladruno_brick]] (`-geom <linear|corot|finite>`), and any
  later Ladruno solid (tets/hexes) — the layer is element-agnostic, keyed only on
  nodal kinematics + ∇N + the core's local f/K.
- **Reference (copy patterns / theory):**
  - `SRC/coordinateTransformation/CorotCrdTransf3d.cpp` — the beam precedent for
    the localize/globalize + geometric-stiffness pattern (adapt to continuum).
  - Belytschko et al. Ch. 4 (§4.6 corotational; §4.5 Voigt; TL/UL) — derivations;
    consistent corotational K_geo: Rankin–Brogan / Felippa–Haugen (EICR).
  - de Souza Neto–Perić–Owen **Ch. 14** — the *spatial* log-strain multiplicative
    framework the seam-3 adaptor implements: §14.3.2 (Eulerian log strain
    εᵉ = ½ln **b**ᵉ), §14.4.2 (backward **exponential-map** integrator), Box 14.3
    (`MATISU`: the exact adaptor algorithm). **F-bar is Ch. 15 §15.1** (eq. 15.5,
    F̄ = (det**F**₀/det**F**)^⅓ **F**); **EAS is Ch. 15 §15.2**.

## How

### Interface sketch

```cpp
class SolidTransformation {
 public:
  // ---- seam 1: kinematics ledger (the part present in every method) ----
  //   fills the strain measure at a GP from reference grads + nodal disp.
  //   v1 small-strain core wants engineering strain; finite wants F.
  virtual void strainAtGP(const double dNdX[][3], const Vector* u8,
                          Vector& eps /*6*/, Matrix* F /*3x3, finite only*/) = 0;

  // ---- seam 2: localize global -> core frame (corot strips rigid R) ----
  virtual void localizeDisp(const Vector& uGlobal, Vector& uLocal) = 0;

  // ---- seam 3 boundary: globalize core f/K back, add geometric stiffness ----
  virtual void globalizeForce (const Vector& fLocal,  Vector& fGlobal) = 0;
  virtual void globalizeStiff (const Matrix& kLocal,  Matrix& kGlobal,
                               const Vector& stressGP /*for K_geo*/) = 0;

  // ---- recorder frame (contract §B.2): current basis for per-step localAxes ----
  virtual bool frameTimeVarying() const = 0;      // corot/finite -> true
  virtual void currentFrame(Vector& vxyz9) const = 0;
};
```

The element loop becomes method-agnostic. **Seam-3 signature (aligned with the
material-wrapper plan, 2026-06-01):** `finite` passes the **full F** (not F_Δ) to a
`FiniteStrainNDMaterial::setTrialF(F)`; the adaptor stores `F_committed`, derives
`F_Δ`, and holds `J=det F`. `getStress()` returns **Cauchy σ**, `getTangent()`
returns the **spatial tangent a** — so the element treats the adaptor as an ordinary
`NDMaterial` and assembles `∫ Bᵀσ dv` spatially.

```
for each GP:
    transform->strainAtGP(dNdX, u8, eps, &F);     // ledger: eps (linear/corot) | full F (finite)
    if (finite) mat->setTrialF(F);                  // seam 3: full F in; adaptor derives F_Δ, J
    else        mat->setTrialStrain(eps);           // linear/corot: small ε in (de-rotated for corot)
    sigma = mat->getStress();  D = mat->getTangent();// finite: σ=Cauchy, D=spatial a (14.99)
    accumulate fLocal, kLocal with the CORE formulation (std/bbar/uri/eas)  // oblivious
transform->globalizeForce(fLocal, fGlobal);         // + nothing (linear) / + K_geo (corot/finite)
transform->globalizeStiff(kLocal, kGlobal, sigma);
```

### `corot` — the corotational method (v2)

1. From current nodal positions, build best-fit element frame; extract rotation
   **R** by polar decomposition of the deformation gradient (Higham iteration on
   a centroidal **F**, or the element best-fit polar of Crisfield/Rankin).
2. `uLocal` = de-rotated deformational displacement (global minus rigid-body part).
3. Core returns `fLocal`, `kLocal` from the **unchanged** small-strain kernel.
4. `fGlobal = Tᵀ fLocal`; `kGlobal = Tᵀ kLocal T + K_geo`, where **K_geo** comes
   from the dependence of **R**/T on **u** (the corotational geometric stiffness).
5. `frameTimeVarying()=true`; `currentFrame()` returns **R**'s columns for the
   recorder's per-step `LOCAL_AXES` ([[ladruno_element_contract]] §B.2).

Validity: large displacement + **large rotation**, small material strain. Reuses
the full small-strain material library and every core formulation.

### `finite` — spatial multiplicative method (v3)

**Spatial, not total-Lagrangian.** The earlier draft mixed a TL story
(Green–Lagrange **E** → PK2 **S**) with a spatial log-strain adaptor; those are
two *different* finite-strain formulations and don't compose. We adopt the
**spatial multiplicative** route (Simo / Miehe 1998 / de Souza Neto Box 14.3),
because it is the one where the small-strain return map is reused *verbatim*
(de Souza Neto eq. 14.89–14.90 has the identical functional format as the
infinitesimal 7.25). Cauchy stress drops straight into the element internal force.

1. Ledger emits the **full** deformation gradient **F** = I + Σ uₐ⊗∇XNₐ at each
   GP (from reference grads + total nodal disp) — no rotation stripped.
2. Hand **F** to the **seam-3 material adaptor** via `setTrialF(F)` (not the raw
   `setTrialStrain`). The adaptor stores `F_committed`, derives the incremental
   **F**_Δ = **F** **F**_committed⁻¹ and `J=det F` itself, then runs Box 14.3
   (`MATISU`): elastic-trial left Cauchy–Green
   **b**ᵉ_trial = **F**_Δ **b**ᵉₙ **F**_Δᵀ → spatial log strain
   εᵉ_trial = ½ln **b**ᵉ_trial → **unchanged small-strain return map** → Kirchhoff
   **τ** → Cauchy **σ** = J⁻¹ **τ**. (Why "unchanged" holds: log strain + backward
   **exponential-map** plastic update + isotropy — the trifecta of Remark 14.5;
   the exponential map is also what keeps plastic incompressibility *exact*,
   Remark 14.4.)
3. Core assembles ∫ **B**ᵀ **σ** dv with the **spatial** (small-strain-shaped) **B**
   on the *current* geometry; tangent = material spatial modulus + **spatial
   geometric (initial-stress) K_geo** from current **σ**. The material modulus the
   adaptor returns already folds the **log-map derivative** (∂εᵉ/∂**b**ᵉ) and the
   **b**ᵉ_trial push-forward (de Souza Neto §14.5).

Validity: large **strain**. Only as good as the materials seam 3 can feed (one
isotropic model — J2 / Drucker–Prager — through the adaptor first).

### The boundary with the material wrapper (seam 3)

**Owned by the parallel material-wrapper work — this doc only pins the contract.**
`finite` does *not* call `NDMaterial::setTrialStrain(ε)`; it calls the adaptor:

```
mat->setTrialF(F);                  // full deformation gradient F (3x3); adaptor stores
                                    //   F_committed, derives F_Δ = F·F_committed⁻¹, holds J=det F
sigma = mat->getStress();           // Cauchy σ = J⁻¹ τ (6, Voigt)  ← element internal force
a     = mat->getTangent();          // spatial modulus a (14.99), incl. log-map derivative
```

`mat` is a `FiniteStrainNDMaterial` (the material-wrapper plan's base; `setTrialF`
lives there, **not** on vanilla `NDMaterial`). The adaptor internally (Box 14.3
`MATISU`, *spatial*): `bᵉ_trial = F_Δ bᵉₙ F_Δᵀ → εᵉ_trial = ½ln bᵉ_trial →
[unchanged small-strain return map] → τ → σ = J⁻¹ τ`. Equivalence is exact for
**isotropic** materials (Simo; Miehe 1998). The element supplies the **full F**;
the adaptor derives F_Δ/J and supplies Cauchy stress + spatial tangent. **Neither
half is "finite strain" alone** — they ship together. Signature locked with the
material-wrapper plan (their D4 `setTrialF`; D2 J-seam): pass **F** so the adaptor
owns `J=det F` for the Kirchhoff↔Cauchy conversion.

> **The hidden state — pin this explicitly.** The committed per-Gauss-point state
> the adaptor stores / `commitState` / `sendSelf` is **b**ᵉₙ (or εᵉₙ) — the
> *elastic left Cauchy–Green tensor*, **not** a plastic multiplier and **not**
> total strain (Box 14.3: **b**ᵉₙ = exp[2εᵉₙ]). The de Souza Neto trick feeds
> εᵉ_trial *as if it were the small-strain "total strain"* into an **unmodified**
> `NDMaterial`, which returns the updated εᵉ. So the adaptor must **own b**ᵉₙ and
> wrap the wrapped material's `commit`/`revert`/`sendSelf`/`recvSelf`. This
> state-ownership coupling — not the stress formula — is the genuinely hard part
> of seam 3.

> **Coordination with the material-wrapper team (their note, folded in).** The one
> place to genuinely *share code* is the low-level **3×3 polar/spectral linear
> algebra** — `corot` needs the rotation **R** (polar decomposition); the adaptor
> needs the eigenvalues λᵢ and eigenvectors (spectral, for ½ln **b**ᵉ and the
> log-map tangent derivative). Coordinate on **one robust shared util** there, not
> at the element-composition level. And **confirm the boundary** (worth a 15-min
> sync): our `corot` is *corotational-for-small-strain* (strips **R**, feeds small
> ε to the unchanged material, adds K_geo) — it is **not** a finite-strain
> kinematic layer. If the material wrapper is itself trying to be a finite-strain
> kinematic layer, we'd be building the same spectral/multiplicative machinery
> twice; the clean split is: `corot` owns rotation-stripping, seam-3 owns the
> multiplicative log-strain promotion.

### Composition with the formulation axis

The transformation wraps an **oblivious** core, so the two axes are *mostly*
orthogonal — clean for `linear`/`corot`, with **two genuine couplings under
`finite`** noted below the table (don't assume them away):

| | `linear` | `corot` | `finite` |
|---|---|---|---|
| `std` | v1 ✓ | v2 | v3 |
| `bbar` | v1 ✓ | v2 | v3 (→ F-bar) |
| `uri` | v1 ✓ | v2 | v3 |
| `eas` | v2 | v2/v3 | v3 (large-strain hourglass caution) |

**Coupling 1 — `bbar`+`finite` = F-bar is *not* a clean wrapper.** It is not "the
bbar kernel run inside an oblivious finite wrapper"; it is a **seam-1
(kinematics-ledger) modification of the deformation gradient fed to the
material**: F̄ = (det**F**₀/det**F**)^⅓ **F** (de Souza Neto §15.1, eq. 15.5 —
centroid dilatation **F**₀, i.e. *mean J* not mean dilatation; this is the
geometrically-nonlinear extension of B-bar, works with *any* material unlike
selective reduced integration). So in `finite`, the bbar path lives in the
ledger, not as a core wrapped obliviously. Flagged for v3.

**Coupling 2 — `eas`+`finite` is stability-fragile.** EAS elements have documented
hourglassing/instability under large compressive strain (de Souza Neto §15.2.5);
the `α`-condensation must be re-derived on current geometry. Treat the v3 `eas`
cell as research, not a drop-in.

## Risks / open questions

> [!question]
> **Polar decomposition robustness.** Which extraction — centroidal-F Higham
> iteration vs. element best-fit polar (Crisfield/Rankin)? Best-fit is more robust
> on distorted/warped hexes but costlier. Recommend the **element-independent
> corotational (EICR)** projector (Rankin–Brogan / Felippa–Haugen): robust *and*
> it hands you the consistent tangent machinery (see next item). Decide in v2.
> **These two choices compose, not compete:** the **iterative (Higham) polar** is
> the R-extraction *primitive* (settled below — it sidesteps eigen-of-**U**, so
> repeated stretches never bite us); EICR is the *framing* that wraps it for the
> consistent tangent. v2 decision is best-fit-**F** scope, not the primitive.

> [!question]
> **Shared 3×3 linear-algebra util (cross-team) — now thinner than first stated.**
> Once `corot` commits to the **iterative Higham polar** for **R** (no eigen-of-**U**,
> settled below), the two efforts no longer share *one* polar/spectral kernel:
> `corot` runs an independent polar iteration, while the seam-3 adaptor needs a
> **degeneracy-safe symmetric 3×3 eigensolver** (orthonormal basis at repeated
> eigenvalues — for ½ln **b**ᵉ and the `DISO2` tangent). The genuine shared assets
> shrink to **(a)** that eigensolver, if `corot`'s best-fit path ends up wanting
> spectral **R** after all, and **(b)** the **degenerate-stretch test fixtures**
> (uniaxial/hydrostatic/axisymmetric) both batteries must run. Confirm in the 15-min
> sync *which* of (a)/(b) is actually co-owned — it may be just the test fixtures.

> [!question]
> **Consistent geometric stiffness for `corot`.** The exact K_geo from ∂R/∂u is
> intricate; an approximate/secant K_geo may hurt Newton convergence. Verify
> quadratic convergence on a large-rotation benchmark before shipping.

- **Explicit + corot:** the per-step frame extraction is in the explicit hot
  loop — keep `currentFrame`/polar cheap; cache where the frame is reused.
- **Recorder per-step frame:** `frameTimeVarying=true` makes the recorder write a
  `LOCAL_AXES` time series — confirm the cost is acceptable for large solid meshes
  (it is per-element-per-step).
- **`finite` material coverage:** worthless without seam-3 materials. Sequence v3
  *after* the material-wrapper work lands at least one isotropic model
  (J2/Drucker-Prager) through the log-strain adaptor.
- **Seam-3 state ownership:** the adaptor owns committed **b**ᵉₙ per GP and must
  wrap the material's `commit`/`revert`/`sendSelf`/`recvSelf` (see *The hidden
  state* above). This — not the stress formula — is the integration risk.
- **Spatial consistent tangent (seam 3):** must include the **log-map derivative**
  ∂εᵉ/∂**b**ᵉ and the **b**ᵉ_trial push-forward (de Souza Neto §14.5); a naïve
  "material tangent only" tangent breaks Newton quadratic convergence under finite
  strain.
- **Repeated principal stretches — ownership split (settled).** The 0/0 in the
  log-map tangent, `(ln bᵢ − ln bⱼ)/(bᵢ − bⱼ)` for `i≠j`, is the off-diagonal term
  in the derivative of the isotropic tensor function (Miehe 1998 `DISO2`; de Souza
  Neto App. A.5) — i.e. **inside the seam-3 adaptor's `MATISU`**. The fix (branch on
  multiplicity all-distinct / two-equal / three-equal, analytic Taylor limit in the
  degenerate branches — preferred over Miehe's ε-perturbation) is **the
  material-wrapper team's**, not ours. Our three boundary obligations: (1) **don't
  re-introduce it in `corot`** — extract **R** with the iterative (Higham) polar,
  never eigen-of-**U**, so repeated stretches are a non-issue on our side (**R** is
  unique for det**F**>0); (2) our `finite` T0/T1 battery is **the net that catches
  their bug** — it must include uniaxial (2-fold), hydrostatic/equibiaxial (3-fold)
  and axisymmetric cases *specifically* to exercise the multiplicity branches, since
  those degeneracies are load-bearing from the first verification run; (3) the
  seam-3 `update()` contract obligates finite, non-NaN **σ** + tangent under
  repeated stretches — their guarantee, our test.

## Roadmap

- **v1 ([[09_ladruno_brick]]):** ship `SolidTransformationLinear` only (identity);
  prove the seam by routing the small-strain brick through it with zero behaviour
  change (regression vs. a direct-kernel build to ~1e-12).
- **v2:** `SolidTransformationCorot` + the rigid-rotation patch test (zero stress
  under arbitrary rotation) + large-rotation Newton-convergence benchmark; wire
  `frameTimeVarying`/`currentFrame` into the recorder.
- **v3:** `SolidTransformationFinite` (**spatial multiplicative**, Box 14.3) against
  the seam-3 log-strain adaptor; exponential-map update + spatial tangent (§14.5);
  F-bar for `bbar` (§15.1, eq. 15.5); necking-of-a-bar (§14.9.2) + geotech
  large-deformation benchmarks; must reduce to `corot` then `linear` as strain→0.

**Phase 0 (gates v3, no element code):** lock the seam-3 signature with the
material-wrapper team — *spatial route* (b ᵉ / Kirchhoff τ / Cauchy σ), the
**state-ownership** contract, and the **shared 3×3 polar/spectral util**. Pin
before either side hardens.

## Implementation log

- 2026-05-31 — Plan drafted. Named *solid transformation* (not "corot") because
  corotational is one of three hosted methods (linear/corot/finite). Established
  the framing: kinematics ledger is **intrinsic** (present in every method); the
  "wrapper" is a discipline for *where* that computation lives so the constitutive
  core stays oblivious. Interface sketch + the three methods + the orthogonal
  composition with the `std/bbar/uri/eas` formulation axis. Pinned the **seam-3
  boundary** with the parallel material-wrapper work (element supplies F; adaptor
  returns S/ℂ via log-strain; exact for isotropic). Solid analogue of beam
  `CorotCrdTransf` (which is frames-only upstream). First consumer
  [[09_ladruno_brick]] ships `linear` only; corot=v2, finite=v3.
- 2026-06-01 — Assessed against the reference library (de Souza Neto–Perić–Owen
  Ch. 14/15, read directly). Corrections folded in: **(1)** `finite` relabeled
  *total-Lagrangian* → **spatial multiplicative** — the draft mixed Green–Lagrange
  **E**/PK2 **S** (TL) with a spatial log-strain adaptor; they don't compose. Adopt
  the spatial route (Box 14.3 `MATISU`): **b**ᵉ_trial=**F**_Δ**b**ᵉₙ**F**_Δᵀ →
  εᵉ=½ln**b**ᵉ → unchanged small-strain return map → τ → σ=J⁻¹τ. **(2)** Named the
  **exponential-map** integrator (eq. 14.72) as the load-bearing ingredient that
  makes the small-strain return map reusable verbatim (Remark 14.5) and plastic
  incompressibility exact (Remark 14.4). **(3)** Pinned seam-3 **hidden state**:
  the adaptor owns committed **b**ᵉₙ per GP and wraps the material's
  commit/revert/sendSelf — the real integration risk, not the stress formula.
  **(4)** Spatial consistent tangent needs the log-map derivative + push-forward
  (§14.5). **(5)** Citation fixes: log-strain plasticity = Ch. 14; **F-bar = Ch. 15
  §15.1** (eq. 15.5, mean *J*); EAS = §15.2. **(6)** Two `finite` couplings that
  break strict axis-orthogonality (bbar+finite=F-bar is a ledger mod; eas+finite is
  stability-fragile §15.2.5). **(7)** Material-team coordination: share **one 3×3
  polar/spectral util** (corot's **R** = adaptor's λᵢ/eigenvectors), not at the
  composition level; **confirm `corot` is corotational-for-small-strain, not a 2nd
  finite-strain kinematic layer** (else both teams build the same machinery — 15-min
  sync). Added a **Phase 0** gate (lock seam-3 signature before v3 code on either
  side).
- 2026-06-01 — Element-collision check **resolved**: `LadrunoBrick` is the sole
  solid element; the material-wrapper plan (`09_finite_strain_material_wrapper.md`)
  explicitly builds *no* element. Seam-3 signature **aligned** to that plan: element
  passes the **full F** to `FiniteStrainNDMaterial::setTrialF(F)` (was F_Δ); adaptor
  derives F_Δ/J and returns Cauchy σ (`getStress`) + spatial tangent **a**
  (`getTangent`). Updated the interface sketch, element loop, `finite` steps 1–2,
  and the seam-3 boundary block to match. v1 code starting now.
- 2026-06-02 — **F-bar (`bbar` + `finite`) IMPLEMENTED + locked against the primary
  source.** Coupling 1 (above) realized: F-bar lives in the element kinematics
  ledger, not the transformation. `SolidTransformationFinite` stays a pure
  identity; `LadrunoBrick` branches on `formulation==BBAR` inside the finite path.
  Locked spec (de Souza Neto, Perić & Owen 2008 Ch.15 §15.1, read directly — the
  book is on `Desktop\FEM expert\books\`):
  - **Modified gradient (eq 15.5):** `F̄ = (det F₀/det F)^(1/3) F`, F₀ = deformation
    gradient at the element centroid (ξ=0). `det F̄ = J₀`, isochoric part = the GP's
    (Remark 15.1). New per-element work: one extra `shp3d` at the centroid for J₀
    (and G₀ for the tangent) — `LadrunoBrick::centroidFbar()`.
  - **Residual (eq 15.9 / Remark 15.2):** UNCHANGED — standard B on the *actual*
    configuration, `dv = J dV` with the **actual** J (not J₀); only the stress input
    changes (material driven by F̄ ⇒ σ̄). So `formResidAndTangentFinite`'s residual
    code needed no edit; `updateFinite` just scales F→F̄ before `setTrialF`.
  - **Tangent (Prop 15.1, eq 15.10/15.11):**
    `K = ∫Gᵀ a|_{F̄} G dv + ∫Gᵀ q (G₀−G) dv`, `q = (1/3) a:(I⊗I) − (2/3)(σ̄⊗I)`,
    with `a` the **full** spatial tangent (= `c − σ̄δ`, the FD-verified modulus the
    std term already uses). In code: `M_ij = (1/3)Σₚ a4_ijpp − (2/3)σ̄_ij`, then
    `K_{(a,i)(b,k)} += (Σⱼ g_{a,j} M_ij)(G₀_{b,k} − g_{b,k}) dv`.
  - **The coefficient was the trap.** A 3-derivation adversarial workflow (blind,
    multi-lens) + my own first pass all produced `q_ij = (1/3) c_ijmm` via the
    spatial shortcut `dσ̄ = c̄:sym(L̄)` — **missing the `−(2/3)σ⊗I`**. The shortcut
    omits the J-factor / reference-config terms the book's full derivation
    (eq 15.12–15.17, via P̂ and the `(detF₀/detF)^(−2/3)` factor) captures. The
    workflow's own FD "confirmation" left a ~0.03 residual it wrote off as a crude-c̄
    artifact — that residual *was* the missing term. Resolved by the primary source;
    the OpenSees FD-tangent test (LogStrain's analytic tangent) is the clean arbiter.
  - **UNSYMMETRIC** in general (book, after eq 15.10) ⇒ `system FullGeneral`; the
    F-bar test does not assert symmetry. Logged in [[LEDGER_quirks]].
  - Parser: `-geom finite -formulation bbar` now accepted (was rejected); uri/eas +
    finite still reserved. Tests: `tests/test_ladrunoBrick_finite.py` adds an
    unsymmetric-aware FD-tangent gate, reduce-to-std-on-homogeneous-F, and the **T4**
    volumetric-locking cantilever (std locks as ν→0.5, F-bar stays compliant).
