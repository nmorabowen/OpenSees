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
work — referenced here as *seam 3*, see *The boundary with the material wrapper*).

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
| **`finite`** (TL) | large **strain** | F → E (Green–Lagrange) | yes (initial-stress) | finite-strain *via* seam-3 log-strain adaptor |

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
  - Belytschko et al. Ch. 4 (§4.6 corotational; §4.5 Voigt; TL/UL) — derivations.
  - de Souza Neto–Perić–Owen Ch. 14 — F-bar & the logarithmic-strain framework
    that the material adaptor implements.

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

The element loop becomes method-agnostic:

```
for each GP:
    transform->strainAtGP(dNdX, u8, eps, &F);     // ledger
    material->setTrialStrain(eps);                  // (or adaptor->update(F) — seam 3)
    sigma = material->getStress();  D = material->getTangent();
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

### `finite` — total-Lagrangian method (v3)

1. Ledger emits the full **F** = I + Σ uₐ⊗∇XNₐ at each GP (no rotation stripped).
2. Hand **F** to the **seam-3 material adaptor** (not the raw `setTrialStrain`):
   it computes E_log, calls the small-strain return map, returns PK2 **S** /
   Kirchhoff **τ** and the consistent tangent.
3. Core assembles with the displacement-dependent **B = B_L + B_NL(u)**; tangent
   = material part + **initial-stress (geometric) K_geo** from current **S**.

Validity: large **strain**. Only as good as the materials seam 3 can feed.

### The boundary with the material wrapper (seam 3)

**Owned by the parallel material-wrapper work — this doc only pins the contract.**
`finite` does *not* call `NDMaterial::setTrialStrain(ε)`; it calls the adaptor:

```
adaptor->update(F);                 // F in (3x3)
S   = adaptor->getPK2();            // or Kirchhoff τ
Dc  = adaptor->getMaterialTangent();// ∂S/∂E (consistent, incl. log-map geometry)
```

The adaptor internally: `F → C=FᵀF → spectral → E_log = ½ln C →
[unchanged small-strain return map] → T → S, ℂ`. The equivalence is exact for
**isotropic** materials (Miehe 1998). The element supplies **F**; the adaptor
supplies stress+tangent. **Neither half is "finite strain" alone** — they ship
together. Lock this signature before either side hardens.

### Composition with the formulation axis

The transformation wraps an **oblivious** core, so the two axes are orthogonal:

| | `linear` | `corot` | `finite` |
|---|---|---|---|
| `std` | v1 ✓ | v2 | v3 |
| `bbar` | v1 ✓ | v2 | v3 (→ F-bar) |
| `uri` | v1 ✓ | v2 | v3 |
| `eas` | v2 | v2/v3 | v3 (large-strain hourglass caution) |

(`bbar` under `finite` becomes **F-bar**, de Souza Neto 1996 — mean *J* instead
of mean dilatation; flagged for v3.)

## Risks / open questions

> [!question]
> **Polar decomposition robustness.** Which extraction — centroidal-F Higham
> iteration vs. element best-fit polar (Crisfield/Rankin)? Best-fit is more robust
> on distorted/warped hexes but costlier. Decide for `corot` in v2.

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
- **Repeated eigenvalues** in the log-strain spectral decomposition (seam 3) — the
  adaptor's problem, but the element must not assume distinct stretches.

## Roadmap

- **v1 ([[09_ladruno_brick]]):** ship `SolidTransformationLinear` only (identity);
  prove the seam by routing the small-strain brick through it with zero behaviour
  change (regression vs. a direct-kernel build to ~1e-12).
- **v2:** `SolidTransformationCorot` + the rigid-rotation patch test (zero stress
  under arbitrary rotation) + large-rotation Newton-convergence benchmark; wire
  `frameTimeVarying`/`currentFrame` into the recorder.
- **v3:** `SolidTransformationFinite` (TL) against the seam-3 log-strain adaptor;
  F-bar for `bbar`; geotech large-deformation benchmark.

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
