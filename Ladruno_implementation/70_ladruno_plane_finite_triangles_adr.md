---
title: "ADR 70 — Plane finite-strain continuum family: shared 2D finite-strain kernel + LadrunoCST -geom finite + new LadrunoLST (T6)"
status: in progress — P0 shipped (#540), P1 shipped (LadrunoQuad -geom finite, #543 + review fixes #545), P2 (LadrunoCST -geom finite) landing; NEXT P3 (LST)
---

# ADR 70 — Plane finite-strain triangles on a shared 2D finite-strain kernel

**Status:** draft. Planning only — no code yet. Extends the plane-elements line
([[25_ladruno_plane_elements_adr]]) and consumes the finite-strain material seam
just shipped as `LogStrain2D` (33016, ADR-25 P5 material half, #536). The frontier
thesis [[26_ladruno_plane_frontier_adr]] frames *why* the triangle is hard; this
ADR is the *engineering* plan for the finite-strain (large-strain) axis on the 2D
continuum family, and the introduction of a proper linear-strain triangle.

Companion theory loaded for this ADR: `/fem-mechanics-expert` (dSNPO 2008 Ch. 14–15,
Belytschko Ch. 4/8), `/abaqus` element technology (CPE3 / CPE6 / CPE6M / hybrid),
`/ls-dyna` (mean-dilatation / selective-reduced volumetric control).

---

## 1. Context & problem

`LogStrain2D` (33016) lifts any small-strain 3D material to a genuine **2D
finite-strain** (Hencky) material behind an order-3 PlaneStrain/PlaneStress face
with a `setTrialF(2×2 F)` seam. But **no 2D element drives it yet** — `LadrunoQuad`
(33007) and `LadrunoCST` (33008) are `-geom linear` only (both carry the explicit
`// … -geom linear only (no theGeom)` guard). So the finite-strain plane capability
is stranded on the material side, exactly as the 3D `LogStrainNDMaterial` was before
`LadrunoBrick -geom finite` (#76) caught up.

The task is to build the **element `-geom finite` path** for the 2D continuum
family. Two forces shape the design:

1. **Reuse across a family, not one element.** The finite-strain *mechanics*
   (compute F, F-bar, `setTrialF`, push-forward gradients, assemble the
   updated-Lagrangian residual and the non-minor-symmetric consistent tangent
   `a = c − σ_il δ_jk`) is **identical** for every 2D continuum element — only the
   shape functions / integration rule differ. The intended consumers are
   `LadrunoQuad` (Q4), `LadrunoCST` (T3), a **new** `LadrunoLST` (T6), and later
   `BezierTri6`. Inlining the mechanics four times would quadruple the hardest,
   most error-prone code. → a **shared, element-agnostic 2D finite-strain kernel**.

2. **The CST is honestly a dead end for the target physics** (the ADR-26 verdict,
   restated): constant strain ⇒ (a) **catastrophic** incompressible/J2 locking
   (one isochoric constraint per element, and **element-level b-bar/F-bar cannot
   save a constant-strain simplex — there is nothing to average**), and (b)
   **mesh-direction-biased** localization (a constant strain field cannot carry a
   shear band at an arbitrary angle). So the fork needs a **linear-strain triangle
   (LST / T6)** as the *usable* triangular substrate for finite-strain plasticity —
   the T6 strain field is linear, so it carries an inclined band and, crucially,
   **has a non-constant J that element-level F-bar CAN average.**

This ADR therefore delivers: (i) the shared kernel; (ii) `-geom finite` on the two
existing continuum elements; (iii) a new `LadrunoLST`; and (iv) a clear, theory-
grounded **per-element F-bar strategy**.

> [!note] Relationship to ADR-25 P4/P5
> ADR-25 scoped `-geom corot` (P4) and `-geom finite` (P5) as future phases of
> `LadrunoQuad`/`LadrunoCST`. This ADR **realizes the `-geom finite` half** for the
> whole family on the shared kernel and adds `LadrunoLST`. **corot (P4) stays out of
> scope** — corot is large-*rotation*/small-*strain* with rigid-rotation stripping,
> a different concern best served by a `SolidTransformation2DCorot` (the 2D analogue
> of the 3D `SolidTransformationCorot` EICR seam), not the finite-strain kernel.

---

## 2. Decision (summary)

1. **Build a shared, element-agnostic 2D finite-strain assembly kernel** in
   `SRC/element/ladrunoPlane/`:
   - `LadrunoFiniteStrain2DKernel.h` — **pure** (plain doubles, no OpenSees deps),
     unit-testable with g++ against a numpy oracle exactly like `LogStrainKernel.h`.
     Holds: F from reference gradients, F-bar scale (power **½**), F⁻¹ push-forward
     to spatial gradients, `a = c − σ_il δ_jk`, the residual and tangent
     contractions, and the F-bar G₀ coupling.
   - `LadrunoFiniteStrain2D.{h,cpp}` — thin OpenSees glue: given per-GP reference
     shape gradients (+ centroid gradients for F-bar), the nodal-disp `Vector`,
     `FiniteStrainND2DMaterial*` pointers, thickness and GP weights, it runs the
     `setTrialF` update loop and fills the `resid` `Vector` / `stiff` `Matrix`.
     Element-agnostic (any node count, any GP rule).

2. **`LadrunoQuad -geom finite`** (Q4) — first consumer / proving ground:
   `-formulation {std|bbar}` × `-geom finite`, with **2D F-bar** (`bbar`+`finite`).
   Completes ADR-25 P5 for the quad. (`ssp`/`eas`+finite stay reserved as in ADR-25.)

3. **`LadrunoCST -geom finite`** (T3) — `std` only. **F-bar is a no-op on T3**
   (constant J), so v1 CST-finite is the honest over-stiff baseline; the real
   volumetric cure for the simplex is **F-bar-Patch** (nodal-averaged dilatation,
   dSNPO §15.1.9), scoped as a later phase (P4) because it is a cross-element
   (nodal-patch) operation the element cannot do alone.

4. **New `LadrunoLST`** (T6, `ELE_TAG_LadrunoLST = 33016`, ELE registry — a
   different registry from `ND_TAG_LogStrain2D = 33016`, no collision):
   6-node quadratic triangle, **3-point interior integration** (matches upstream
   `SixNodeTri`), `-formulation {std|bbar}` × `-geom {linear|finite}`. Reduces to
   upstream `SixNodeTri` at `-geom linear -formulation std` (the reduce-to gate,
   mirroring Quad↔FourNodeQuad / CST↔Tri31). **Element-level F-bar works** on T6
   (non-constant J) — this is the element that makes triangular finite-strain
   plasticity actually usable.

5. **F-bar strategy per element** (the load-bearing formulation decision):
   | Element | strain field | element-level F-bar | volumetric cure |
   |---|---|---|---|
   | Q4 `LadrunoQuad` | bilinear (J varies) | **works** | `bbar`+`finite` (this ADR) |
   | T6 `LadrunoLST`  | linear (J varies)  | **works** | `bbar`+`finite` (this ADR) |
   | T3 `LadrunoCST`  | constant (J const) | **no-op**  | F-bar-Patch, nodal (P4) |

---

## 3. Formulation (grounded in theory)

### 3.1 Updated-Lagrangian 2D finite-strain assembly (the shared core)

Per Gauss point the element supplies the **reference** shape gradients
`∂Nₐ/∂X` (a = 1…nNodes). The kernel computes, in 2D (indices i,j,k,l ∈ {1,2}):

- **Deformation gradient** (direct / index / Voigt-free — it is a full 2-tensor):

  **F** = **I** + Σₐ **uₐ** ⊗ ∂Nₐ/∂**X**  ⟺  F_iJ = δ_iJ + Σₐ (uₐ)_i ∂Nₐ/∂X_J,   J = det **F** (2×2).

- **Feed the material** the total 2×2 F: `fsm->setTrialF(F)` →
  Cauchy **σ** (in-plane, order-3 Voigt {11,22,12}) + the full in-plane spatial
  modulus **c** via `getSpatialTangentTensor2D(c[2][2][2][2])`.

- **Spatial gradients** (push-forward to the current config): ∂Nₐ/∂x_j = Σ_k ∂Nₐ/∂X_k (F⁻¹)_kj.

- **Internal force** (current config, thickness t, area weight w·detJ₀):

  f_{a,i} = ∫ σ_ij ∂Nₐ/∂x_j t dA  ≈  Σ_gp (σ_ij g_j^a) · J · t · w · detJ₀.

- **Consistent tangent** — the **locked seam-3 contract** (dSNPO eq. 14.99):

  a_ijkl = c_ijkl − σ_il δ_jk,   K_{(a,i)(b,k)} = ∫ (∂Nₐ/∂x_j) a_ijkl (∂N_b/∂x_l) t dA.

  The geometric/initial-stress term −σ_il δ_jk is **NOT minor-symmetric in (k,l)**,
  so it must enter the *full-gradient* contraction (not the conventional ∫ Gᵀ Σ G
  alone). This is the same non-symmetric-c cancellation that makes **K** symmetric,
  already proven for the 3D `LadrunoBrick` finite path — the 2D kernel is a faithful
  i,j,k,l ∈ {1,2} restriction of `LadrunoBrick::formResidAndTangentFinite`.

### 3.2 F-bar and the simplex problem (why T3 ≠ T6)

Near-incompressible / J2-plastic flow imposes J ≈ const. Displacement elements
over-constrain this (volumetric locking). The F-bar cure (dSNPO §15.1, Abaqus
mean-dilatation / B-bar §3.2.4) replaces the dilatation with an **element-average**:

  **F̄** = (J₀ / J)^{1/n} **F**,   with n = **2** in plane (NOT the brick's ⅓ — a copy-paste of the ³√ is a bug),

so every GP shares the centroid dilatation J₀ = det **F**(centroid). The residual is
unchanged except through σ̄(F̄); the tangent gains a generally **unsymmetric** G₀
coupling to the centroid gradient operator (dSNPO eq. 15.10) — the kernel carries it.

- **Q4 / T6**: J varies over the element ⇒ (J₀/J) ≠ 1 ⇒ F-bar is **effective**.
- **T3 (CST)**: J is **constant** over the element ⇒ J₀ ≡ J ⇒ **F̄ = F**, a **no-op**.
  A constant-strain simplex has "nothing to average" (ADR-26; Belytschko Ch. 8).
  The correct cure is **F-bar-Patch** (dSNPO §15.1.9): compute a **nodal** patch-
  averaged J̄ by volume-averaging J over the elements sharing each node, then
  F̄ = (J̄/J)^{1/2} F per element. This is a *nodal/cross-element* operation → it
  cannot live in a single element's Gauss loop; it needs a domain-level J̄ assembly
  pass. Deferred to **P4** (see §6), keeping CST-finite v1 an honest baseline.

### 3.3 The linear-strain triangle (LST / T6)

Six-node quadratic triangle: corner + midside nodes, quadratic N, **linear** strain.

- **Shape functions / integration:** quadratic isoparametric N; **3-point interior
  rule** (barycentric (⅔,⅙),(⅙,⅔),(⅙,⅙), w = ⅙) — integrates the quadratic strain
  energy exactly and matches upstream `SixNodeTri` (Kokot 2020). (The alternative
  3-midside rule is less accurate — `SixNodeTri` ships the interior rule; we match.)
- **Why it fixes CST:** the linear strain field (a) carries an inclined shear band
  (no mesh-line snapping) and (b) has a non-constant J, so **element-level F-bar
  works** — T6 is the F-bar-friendly triangle the family was missing.
- **Residual T6 caveat (from Abaqus):** the plain quadratic triangle (CPE6) still
  locks somewhat near full incompressibility and has the uniform-pressure corner-
  load anomaly; Abaqus ships the *modified* CPE6M (+ hybrid CPE6MH) for this. Our
  answer is the same in spirit: **`bbar`+`finite` F-bar** on T6 for the volumetric
  path; a hybrid u/p T6 is a frontier item, not this ADR.

### 3.4 Cross-code positioning

| Concern | Abaqus | LS-DYNA | This ADR |
|---|---|---|---|
| Linear tri | CPE3/CPS3 (over-stiff, "non-critical regions") | const-strain tri | `LadrunoCST` (honest baseline) |
| Quadratic tri | CPE6 / CPE6M(H) modified+hybrid | higher-order tri | `LadrunoLST` (T6) + F-bar |
| Volumetric cure | B-bar / mean-dilatation §3.2.4; hybrid u/p | selective-reduced / mean-dilatation | F-bar (Q4/T6), F-bar-Patch (T3, P4) |
| Finite-strain update | D-integrated, Hughes-Winget objective rate | Jaumann rate | **F-based** via `FiniteStrainND2DMaterial` (spatial multiplicative Hencky) — *diverges* (total F, not rate) |

Note the deliberate divergence: we drive finite strain by the **total F** through the
Hencky material seam (dSNPO Box 14.3), not by an objective **rate** update of Cauchy
stress (Abaqus/LS-DYNA). This reuses the entire small-strain material library
verbatim and is exact under superposed rotation by construction.

---

## 4. Architecture

```
SRC/element/ladrunoPlane/
  LadrunoFiniteStrain2DKernel.h     # PURE doubles; g++-testable vs numpy oracle
  LadrunoFiniteStrain2D.{h,cpp}     # OpenSees glue: update() + resid/stiff fill
  LadrunoQuad.{h,cpp}               # + Geom axis, finite path -> kernel   (33007)
  LadrunoCST.{h,cpp}                # + Geom axis, finite path -> kernel   (33008)
  LadrunoLST.{h,cpp} + OPS_*.cpp    # NEW T6 element                       (33016)
```

**Element responsibilities (thin, per element):** shape functions + reference
gradients per GP + centroid; GP rule + weights; thickness; hold materials as
`FiniteStrainND2DMaterial*`; parser `-geom finite`; route `update()` /
`getResistingForce()` / `getTangentStiff()` to the kernel when `isFinite()`.

**Kernel responsibilities (shared, written once):** everything in §3.1–3.2 — F,
F-bar (½), F⁻¹ push-forward, `a = c − σδ`, residual + tangent contraction, G₀ term.

**Why a kernel, not the 3D `SolidTransformation` seam:** the 3D `SolidTransformation`
is a *per-element-instance virtual wrapper* whose real job is corot rotation-
stripping; for pure finite strain it degenerates to an identity globalize while the
*element* does the F/K_geo work. A shared **kernel** (the fork's `*Kernel.h` idiom:
`LogStrainKernel`, `LadrunoJ2Kernel`, `LadrunoEmbeddedKernel`) gives the cross-element
reuse the finite path actually needs, and stays g++/oracle-testable.

**Future consumers:** `BezierTri6` (33000) adopts the kernel by `#include`-ing it and
supplying its Bézier gradients — **zero new mechanics** (ties into
[[project_bezier_transformation_connect]]).

---

## 5. Class tags & registration

- `ELE_TAG_LadrunoLST = 33016` (ELE registry; free per the ADR-66 note "33016-33019
  remain free ELE slots"; per-registry, so no collision with `ND_TAG_LogStrain2D`).
- Kernel + `LadrunoFiniteStrain2D` glue: **no classTag** (plain helper, like
  `LogStrainKernel` / `LadrunoEmbeddedKernel`).
- Vanilla touches (additive, `// Ladruno`): `SRC/classTags.h` (LST tag),
  `FEM_ObjectBrokerAllClasses.cpp` (LST broker), the element command registries
  (Tcl + Py `element LadrunoLST …`), `SRC/element/ladrunoPlane/CMakeLists.txt`.
- `LadrunoQuad`/`LadrunoCST` `sendSelf`/`recvSelf` gain the `Geom` enum int
  (default 0 = linear ⇒ byte-identical for existing models).

---

## 6. Phases & exit gates

| Phase | Scope | Gate |
|---|---|---|
| **P0** | shared kernel `LadrunoFiniteStrain2DKernel.h` + numpy oracle | oracle reproduces: homogeneous-F patch σ/tangent; FD consistent-tangent (unsym-aware) to 1e-7; reduce-to-linear at ε→0; F-bar ½-power vs a near-incompressible closed form; C++-kernel-vs-oracle to 1e-9 (g++, no OpenSees build) |
| **P1** | `LadrunoQuad std+bbar -geom finite` (ADR-25 P5 for the quad) | element FD consistent-tangent (unsym-aware); homogeneous-F patch vs oracle; small-strain reduce-to-`-geom linear`; det F≤0 step-cut; **volumetric-locking cantilever** — `std` locks at ν→0.5, `bbar`(F-bar) does not; end-to-end `LogStrain2D` over `LadrunoJ2`/elastic |
| **P2** | `LadrunoCST std -geom finite` | homogeneous-F patch vs oracle; reduce-to-linear; det-F step-cut; documents (and pins) the F-bar no-op + over-stiff volumetric response as the honest baseline |
| **P3** | new `LadrunoLST` (T6): `std+bbar` × `-geom {linear,finite}` | **reduce-to `SixNodeTri`** at `linear/std` (~1e-9); rank/3-rbm; inclined-band (localizes off mesh lines, vs CST); finite patch vs oracle; **F-bar effective** on T6 (near-incompressible cantilever: `std` locks, `bbar` doesn't) |
| **P4** *(optional / frontier)* | **F-bar-Patch** nodal dilatation for simplices (CST + future tets) | mesh-objective near-incompressible T3 patch (dSNPO §15.1.9); domain-level J̄ pass; compare CST-patch vs LST — quantify whether the patch closes the gap enough to matter |

Adversarial-gate policy (per [[feedback_adversarial_gate_when]]): **full Opus gate at
P0** (novel shared math: the unsym F-bar G₀ term + the non-minor-symmetric tangent
contraction) and **P3** (new element + reduce-to-upstream correctness); tests carry
P1/P2 (mechanical mirrors of the P0 kernel + the shipped 3D finite path).

Each phase is one PR off `ladruno`. P1 also closes ADR-25 P5 for the quad.

---

## 7. Risks / open questions

> [!question]
> **Kernel home / Bézier reach.** `ladrunoPlane/` is the natural home, and
> `bezierTriangle/` can `#include` across dirs. If cross-dir include friction shows
> up, promote the pure kernel to a neutral `SRC/element/` location. Decide when the
> first non-`ladrunoPlane` consumer (BezierTri6) actually wires in.

> [!question]
> **F-bar-Patch appetite (P4).** The nodal-patch dilatation is the *only* thing that
> makes a plain T3 usable near incompressibility, but it is genuinely heavier
> (domain-level J̄ assembly, patch bookkeeping, an added state pass). Given `LadrunoLST`
> already delivers a usable F-bar triangle in P3, is CST-patch worth building, or does
> CST stay the deliberate baseline and LST become the recommended triangle? Lean:
> **ship P0–P3, decide P4 by demand** (matches the ADR-26 "CST is a dead end" stance).

> [!question]
> **LST corot.** This ADR does the *finite* (large-strain) axis. Large-rotation/
> small-strain corot for T6 (midside-node polar/EICR) is deferred with the ADR-25 P4
> `SolidTransformation2DCorot`. Confirm corot stays out of scope here.

- **F-bar power ½, not ⅓** (plane dilatation J = det 2×2 F). Pin with a
  near-incompressible plane-strain cantilever; copying the brick's ^{1/3} is a bug.
- **T6 distortion fragility** (Abaqus caveat): midside-node placement and heavy
  distortion degrade the quadratic triangle; document, and lean on straight-edged
  meshes for the gates.
- **ndf consistency:** continuum nodes are ndf = 2; finite adds no drilling DOF.
- **Backwards compatibility:** new element (LST) + new `-geom` enum defaulting to
  linear ⇒ zero impact on existing models; `sendSelf` gains one int.

---

## 8. References

- de Souza Neto, Perić & Owen (2008), *Computational Methods for Plasticity*:
  Box 14.3 (MATISU spatial multiplicative), §14.5 (14.99 spatial tangent),
  **§15.1 F-bar**, **§15.1.9 F-bar-Patch (simplex)**, §15.2 EAS, §15.3 mixed u/p.
- Belytschko, Liu, Moran & Elkhodary (2014), *Nonlinear FE for Continua & Structures*:
  Ch. 4 (updated Lagrangian, Voigt §4.5), Ch. 8 (locking, B-bar, patch test).
- Abaqus Theory Guide §3.2 (continuum): §3.2.4 B-bar/mean-dilatation, §3.2.6 triangles/tets;
  CPE3/CPS3, CPE6/CPS6, CPE6M(H) modified/hybrid triangles.
- LS-DYNA Theory: selective-reduced / mean-dilatation volumetric control; objective-rate finite-strain update (contrast to our F-based Hencky path).
- Fork: [[25_ladruno_plane_elements_adr]] (plane family), [[26_ladruno_plane_frontier_adr]]
  (CST verdict + triangle frontier), [[09_finite_strain_material_wrapper]] +
  [[project_finite_strain_wrapper]] (`LogStrain2D` / `LogStrainNDMaterial`),
  `LadrunoBrick -geom finite` (the 3D assembly template), upstream `SixNodeTri`
  (T6 reduce-to gate).

## 9. Implementation log

*(filled in as phases land; move to `Ladruno_internal/` when complete)*

### P0 — shared kernel + oracle (#540, merged)

`LadrunoFiniteStrain2DKernel.h` + numpy oracle + C++-vs-oracle diff. See the PR
and [[LEDGER_implementations]]. Reserves `ELE_TAG_LadrunoLST = 33016`.

### P1 — `LadrunoQuad std+bbar -geom finite` (shipped #543, review fixes #545)

Wires the quad (33007) as the kernel's first consumer, closing ADR-25 P5 for the
quad. No new class tag; no vanilla files touched (all edits are fork-owned
`ladrunoPlane/` + the fork test dir).

**What landed**

- **Geometry axis.** New `enum class Geom { LINEAR = 0, FINITE = 2 }` (values mirror
  `LadrunoBrick`), a `geom` member, `-geom linear|finite` parser flag, and `Geom`
  carried in `sendSelf`/`recvSelf` (data vector 14 → 15; default LINEAR ⇒ existing
  models unchanged in value). `Print` reports the axis.
- **`updateFinite()`.** Per GP: reference gradients ∂Nₐ/∂X straight from the existing
  `shapeFunction` (it already differentiates w.r.t. the *reference* nodal coords),
  `F = I + Σ uₐ⊗∂Nₐ/∂X` via the kernel, F-bar scale `(J₀/J)^{1/2}` for `bbar`
  (kernel n = 2, **not** the brick's ⅓), det-F ≤ 0 / J₀ ≤ 0 step-cut guards, then
  `setTrialF(2×2 F)` on the material cast to `FiniteStrainND2DMaterial`.
- **`formFinite(tang_flag)`.** Recomputes the *unbarred* F per GP for the spatial
  gradients `∂Nₐ/∂x = ∂Nₐ/∂X·F⁻¹` and current volume `dv = J·detJ₀·t·w`; assembles
  `f = ∫σ̄ g dv` (kernel `addInternalForce2D`) and, for the tangent, `a = c − σδ`
  (`spatialTangent2D`) → `addTangent2D` + the unsymmetric F-bar G₀ coupling
  (`addFbarCoupling2D`, using the centroid spatial-gradient operator g₀). σ̄ and c
  come from the material set in `updateFinite` (dSNPO Remark 15.2: the F-bar residual
  is the standard resid at σ̄ with the unbarred config).
- **Row-major → column-major.** The kernel writes a row-major `Kloc[64]`; it is copied
  element-wise into the shared column-major `Matrix K` (a blind memcpy would
  **transpose** the genuinely-unsymmetric F-bar tangent — deliberately avoided).
- **`getInitialStiff` unchanged.** It keeps the symmetric reference-config `B̄ᵀD₀B̄`
  seed — a well-posed symmetric operator for Rayleigh βK0 / `-initial` / eigen
  seeding. This is a *deliberate seed*, **not** the true consistent F-bar tangent:
  that stays unsymmetric even at F = I (the F-bar coupling does **not** vanish —
  `q = ½ c:I ≠ 0`, and `g₀ ≠ g` for a distorted map), same rationale `LadrunoBrick`
  documents. (The "coupling vanishes when g₀ = g" wording in the original P1 comment
  was wrong; corrected in #545.) `getMass` (reference-config consistent) also unchanged.
- **Guards.** ssp/eas + finite and **PlaneStress + finite** refused at parse (the
  finite volume weight `dv` omits the thickness stretch λ ⇒ plane-strain only until a
  later phase carries λ + its linearization); a non-`FiniteStrainND2DMaterial` under
  `-geom finite` is rejected at the parser (defensive `dynamic_cast` in
  `updateFinite`/`formFinite` as a backstop); explicit bulk-viscosity is stripped
  (diagnostic) under finite (its force block lives only in the small-strain path).

**Verification.** The finite-strain *assembly math* is already gated building-free by
P0 (the C++ kernel diff covers exactly the Q4 std+F-bar f/K this element assembles).
P1 adds `tests/test_ladrunoquad_finite.py` (openseespy) for the element integration
the build alone can prove: reduce-to-linear (std+bbar collapse to `-geom linear` at
ε→0), a homogeneous finite-stretch patch whose GP Cauchy stress matches the LogStrain2D
oracle σ(F) at the solved F, large-strain Newton convergence over elastic and a 3D J2
inner, F-bar volumetric-locking relief at ν→0.5, and the parse/material guards.

**Post-merge review fixes (#545).** An adversarial pass CONFIRMED the plane-strain
math sound (F-bar ½ + n=2 G₀, `a = c − σδ`, unbarred-J `dv`, the no-transpose copy)
and applied three fixes: (i) **PlaneStress + `-geom finite` refused at parse** — the
finite `dv` omits the thickness stretch λ, so plane stress would carry an ~ν·strain
force error; lifting it needs λ in `dv` + its linearization (a real later phase);
(ii) the shared static `K`/`P` are zeroed on the `formFinite` spatial-tangent error
path (no stale sibling-matrix leak); (iii) the `getInitialStiff` rationale comment
corrected (above). Battery 9 → 11: the added discriminating gates are an
**external-equilibrium assert** (σ × deformed area × t = applied load) and a **Newton
iteration-count ceiling** — a residual-only (converged-result) test cannot catch a
transposed K or a wrong-sign G₀ because Newton still converges to the right answer.

**Next:** P2 = `LadrunoCST std -geom finite` (F-bar is a no-op on the constant-strain
T3 — honest over-stiff baseline); P3 = new `LadrunoLST` (T6).

### P2 — `LadrunoCST std -geom finite`

Wires the T3 (33008) as the kernel's second consumer. No new class tag; no
vanilla files touched (fork-owned `ladrunoPlane/` + the fork test dir + the
banner regen). A mechanical mirror of P1 restricted to the single centroid GP,
per the ADR gate policy (tests carry P2; no separate adversarial gate).

**What landed**

- **Geometry axis.** Same `enum class Geom { LINEAR = 0, FINITE = 2 }` as the
  quad, `-geom linear|finite` parser flag, `Geom` carried in `sendSelf`/`recvSelf`
  (data vector 13 → 14; default LINEAR ⇒ existing models unchanged in value).
  `Print` reports the axis.
- **`updateFinite()` / `formFinite()`.** The T3 is one centroid GP with constant
  reference gradients: F once per element via the kernel, `setTrialF(2×2 F)` on
  the material cast to `FiniteStrainND2DMaterial` (dynamic_cast, graceful fail),
  then the kernel residual `∫σ g dv` + consistent tangent `a = c − σδ` with
  `dv = J·detJ₀·t·w`. **No F-bar lane at all** — on the constant-strain T3,
  J ≡ J₀ so F̄ = F identically (§3.2); CST-finite is the deliberate, honest,
  volumetrically over-stiff baseline. det F ≤ 0 cuts the step through the
  LogStrain2D `setTrialF` guard. Row-major kernel K copied element-wise into the
  column-major shared Matrix (convention kept exact even though the T3 std
  tangent is symmetric); shared static K/P zeroed on the error path (P1 review
  lesson).
- **Guards.** PlaneStress + finite refused at parse (thickness-stretch λ omitted
  from the volume weight — plane-strain only, mirrors the quad); a
  non-`FiniteStrainND2DMaterial` refused at parse + re-checked in `setDomain`
  for broker-built elements; explicit bulk-viscosity stripped (diagnostic) under
  finite. `getInitialStiff` stays the symmetric reference `BᵀD₀B` seed.

**Verification.** `tests/test_ladrunocst_finite.py` (openseespy): reduce-to-linear
at ε→0; homogeneous finite-stretch patch on a two-triangle unit square (both
centroid stresses vs the LogStrain2D oracle σ(F) at the solved F, plus the
deformed-config equilibrium check σxx·(1+b)·t = applied Fx); large-strain Newton
convergence over elastic and a 3D J2 inner; an element-inverting one-step gate
(det F ≤ 0 → analyze fails gracefully, interpreter survives); the **honest-
baseline pin** — the CST/F-bar-quad response ratio collapses as ν→0.5 (locking),
so a future fake per-element T3 "F-bar" (or a quad F-bar regression) trips it;
and the parse guards.

**Next:** P3 = new `LadrunoLST` (T6, 33016) — the F-bar-friendly triangle.
