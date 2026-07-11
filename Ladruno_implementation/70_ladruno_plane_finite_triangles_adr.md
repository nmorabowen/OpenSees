---
title: "ADR 70 — Plane finite-strain continuum family: shared 2D finite-strain kernel + LadrunoCST -geom finite + new LadrunoLST (T6)"
status: in progress — P0 shipped (#540), P1 shipped (LadrunoQuad -geom finite, #543 + #545), P2 shipped (LadrunoCST -geom finite, #550), P3 shipped (LadrunoLST T6 std-only, #554); P4 design spike DECIDED 2026-07-11 (disjoint 2-triangle patch macro-element, exact patch-local unsymmetric tangent — see §9 P4 spike), implementation pending
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
   volumetric cure for the simplex is **F-bar-Patch** (disjoint 2-triangle
   patches sharing J̄ = v_patch/V_patch, dSNPO §15.1.9), scoped as a later phase
   (P4) because it is a cross-element operation a lone element cannot do —
   though the P4 spike showed it fits a patch **macro-element** exactly (§9).

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
   | T6 `LadrunoLST`  | linear (J varies)  | ~~works~~ **REFUTED at P3** | pair macro-element, patch-P1 dilatation (P4b, by demand) |
   | T3 `LadrunoCST`  | constant (J const) | **no-op**  | pair macro-element, patch-constant J̄ per dSNPO §15.1.9 (P4a) |

   > [!warning] P3 amendment (2026-07-10) — the T6 row above was WRONG.
   > "J varies so F-bar can average" is true but insufficient: constant
   > element-mean dilatation (B-bar, or the centroid F-bar scale (J₀/J)^½) is
   > **rank-deficient on the T6**. The two quadratic **conformal** modes
   > (Re/Im z²) have identically zero deviatoric strain, and their linear
   > dilatation has zero element mean — the averaged operator assigns them
   > **zero energy**. A free element shows **5** zero-energy modes (3 RBM + 2
   > spurious; stacked-B̄ rank 7 of the required 9). Verified three ways at P3:
   > analytically, by numpy rank, and by `eigen` on the compiled element. The
   > Q4 is immune only because z²-type fields are not in the bilinear space.
   > Consequence: `LadrunoLST` ships **std-only** (`-formulation bbar`
   > parser-refused with this rationale); the triangle volumetric cure is P4.
   > **P4 spike (2026-07-11) decided the cure** — disjoint 2-triangle patch
   > macro-element (patch-constant J̄ for T3 per dSNPO §15.1.9; patch-P1
   > projected dilatation for T6), exact patch-local unsymmetric tangent.
   > Both candidates named above were refuted as stated: patch-constant J̄
   > does NOT cure the T6 (the conformal modes re-center at the patch
   > centroid), and an element-local P1 projection at the 3-pt rule is the
   > identity (a formulation no-op). See the §9 P4 design-spike log.

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
  The correct cure is **F-bar-Patch** (dSNPO §15.1.9): split the mesh into
  **non-overlapping patches** of simplices (2 triangles suffice in 2D, eq. 15.34 +
  p. 667 guidance) and give every element of a patch the shared patch dilatation
  J̄ = v_patch/V_patch, i.e. F̄ = (J̄/J)^{1/2} F per element. *(Earlier drafts of
  this ADR glossed §15.1.9 as a "nodal" average — wrong: the book method is
  disjoint-patch, and that distinction is exactly what makes its consistent
  tangent patch-local, eqs. 15.37–15.38.)* This is still a cross-element
  operation → it cannot live in a single element's Gauss loop, but because the
  patches are disjoint it CAN live in a patch **macro-element** (P4 decision, §9).
  Deferred to **P4** (see §6), keeping CST-finite v1 an honest baseline.

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
| **P4a** | **T3 pair macro-element** — dSNPO §15.1.9 F-bar-Patch on disjoint 2-triangle patches, exact patch-local unsymmetric tangent (design DECIDED, see §9 spike) | free-pair rank = exactly 3 RBM; FD consistent-tangent on the 8-DOF pair (unsym-aware, at finite stress — the cross blocks vanish at F=I); mesh-objective near-incompressible T3 patch vs the F-bar quad (dSNPO Fig 15.8 localization + necking anchors); conforming-mesh assert; honest-baseline pin flips (pair no longer locks where std T3 does) |
| **P4b** *(optional / by demand)* | **T6 pair macro-element** — patch-P1 projected dilatation over the 2-element patch (novel/off-label; rank + inf-sup evidence in §9 spike) | free-pair rank = exactly 3 RBM incl. distorted + conforming-curved pairs; pressure/checkerboard gate on structured AND alternating-diagonal/unstructured meshes (inf-sup only measured on straight structured so far); FD consistent tangent at finite stress; full adversarial gate (novel math) |

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

> [!question]- **F-bar-Patch appetite (P4). RESOLVED by the 2026-07-11 spike.**
> The premise shifted twice: (a) P3 refuted element-level F-bar on the T6, so the
> patch is no longer just a T3 nicety — it is the family's ONLY triangle
> volumetric cure; (b) the spike showed the cure needs **no** domain-level J̄
> assembly — disjoint 2-triangle patches make the exact tangent patch-local, so
> the whole thing is an ordinary 4-node (T3) / 9-node (T6) macro-element.
> Decision: **build P4a (T3 pair, published method); P4b (T6 pair, novel
> patch-P1) stays by-demand.** See the §9 spike log for the refutations,
> rank/inf-sup evidence, and the rejected tangent strategies.

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

### P3 — new `LadrunoLST` (T6, `ELE_TAG_LadrunoLST` = 33016) + the F-bar-on-T6 refutation

Quadratic N ⇒ **linear** strain ⇒ an inclined band is representable and
convergence is quadratic. Third consumer of the shared kernel. **The headline
P3 finding, however, is negative:** the planned `bbar`/F-bar lane was
**refuted by the rank gate** (see the §2.5 amendment) — `LadrunoLST` ships
**std-only**, and the triangle volumetric cure moves squarely to P4.

**The refutation (how it was caught).** The canonical single-free-element
zero-energy gate (`assert_zero_energy`, locked T0) returned **5** zero modes
for the as-built B-bar T6 (3 RBM + 2 spurious). Diagnosis: the two quadratic
conformal modes (Re/Im z²) satisfy ε_dev ≡ 0 pointwise while their linear
dilatation has zero element mean, so ANY constant-mean averaging (B-bar, or
the centroid F-bar scale) sees zero strain — stacked-B̄ rank is 7 of the 9
required. Confirmed independently by numpy rank and analytically. Per fork
policy (no uncontrolled mechanism ships; cf. SSP's `Kstab`, the ADR-20
refusal), the lane was removed rather than papered over.

**What landed**

- **`LadrunoLST.{h,cpp}` + `OPS_LadrunoLST.cpp`** in `ladrunoPlane/`:
  **std** × `-geom {linear|finite}` × `-type` (finite PlaneStrain-only,
  mirroring the family), `-thick/-rho/-body/-pressure/-bv`. Node order and the
  **3-point interior rule** match upstream `SixNodeTri` (corners n1-n3 CCW,
  midsides n4=(1-2), n5=(2-3), n6=(3-1)); `shapeFunction` is a faithful mirror
  of `SixNodeTri::shapeFunction` (the reduce-to anchor). `-formulation bbar`
  is **parser-refused with the rank rationale** (the message cites the two
  conformal modes and points at P4).
- `getInitialStiff` uses the clean indexed assembly — upstream `SixNodeTri`'s
  version carries FourNodeQuad's 8-stride `matrixData` indexing (latent
  upstream bug, not ported). Consistent quadratic-edge `-pressure` (corner ⅓ /
  midside ⅔ halves) — fixes upstream's side-61 `x4-x6` typo (documented
  divergence). **HRZ lumped mass** (∫ρN² rescaled to the exact total, strictly
  positive) — a second deliberate divergence: upstream's plain N-lumping gives
  exactly ZERO corner masses on the T6 (negative when distorted), unusable for
  explicit dynamics; same decision as the H20 (ADR-72). `charLength = √(2A)`.
- **Finite path (std)**: `updateFinite`/`formFinite` on the kernel — per-GP F,
  `setTrialF`, spatial residual + consistent `a = c − σδ` tangent, det-F
  step-cut via the LogStrain2D guard, error-path K/P zeroing, element-wise
  row-major→column-major copy — the P1/P2 conventions kept (structure mirrors
  `LadrunoCST` with 6 nodes / 3 GPs).
- **Registration** (vanilla, additive `// Ladruno`): `classTags.h` 33016 define
  (replaces the P0 reservation comment; ELE registry, no collision with
  `ND_TAG_LogStrain2D`), broker case, `OpenSeesElementCommands` map +
  classic-Tcl dispatch, `ladrunoPlane/CMakeLists.txt`.
- Note for openseespy cross-checks: upstream `SixNodeTri` registers as
  **`tri6n`** in the shared functionMap (not "SixNodeTri").

**Verification** (`tests/test_ladrunolst_element.py` + `test_ladrunolst_finite.py`;
kernel T6 std math was already gated building-free by the P0 oracle's T6 lane):
reduce-to `SixNodeTri` (~1e-9, PlaneStrain + PlaneStress); rank/3-RBM (std);
**the B-bar refutation pin** (bbar parser-refused — trips if re-enabled without
a rank-sufficient projection); constant-stress patch + **linear-strain-field
reproduction** under a quadratic displacement patch (the CST-impossible field —
the mechanism behind off-mesh-line bands; both on an exact Transformation-
condensed prescribed-field rig with a dummy free DOF — penalty sp constraints
plateau at ~k/α ≈ 1e-7 and cannot meet the locked 1e-9 patch rtol);
reduce-to-linear at ε→0; homogeneous finite-stretch patch vs the LogStrain2D
oracle on a 2-LST square with consistent quadratic-edge loads (⅙, ⅔, ⅙) +
deformed-config equilibrium; finite Newton convergence (elastic + J2);
std-finite Newton iteration ceiling; det-F step-cut; parse guards
(PlaneStress+finite, non-finite material, bbar, bbar+finite, ssp/eas).

**P4 (now load-bearing):** the family's triangle volumetric cure — design
decided by the spike below. Until it ships: quads (`bbar`) where
near-incompressibility matters; triangles are std.

### P4 design spike (2026-07-11) — triangle volumetric cure DECIDED, no code

Design-first spike, as flagged at P3: decide the **tangent strategy** for the
triangle volumetric cure before writing code. Method: dSNPO §15.1.9 read
against the primary text (pp. 665–669), small-strain (F = I) rank analysis in
numpy (`adr70_p4_spike/p4_rank_spike.py`), and a full adversarial gate (novel
math) — an independent twin rebuild of all machinery + distorted/curved/
inf-sup attacks (`p4_adversarial.py`, `p4_infsup_clean.py`,
`p4_curved_conforming.py`). All five spike claims **CONFIRMED** by the gate;
the gate's two apparent refutations were bugs in its own extensions
(a cracked-mesh dof map; a pinv-polluted free-mesh inf-sup), both self-caught.

**DECISION — the tangent strategy.** The cure is built on **disjoint
2-triangle patches implemented as a single macro-element** (union node set:
4 nodes/8 DOF for a T3 pair, 9 nodes/18 DOF for a T6 pair) carrying the
**exact, patch-local, generally unsymmetric consistent tangent** through the
standard element-local assembly interface. dSNPO eqs. 15.37–15.38 give
K^(ee) plus cross blocks K^(es) **only for s in the same patch** — rows/cols
never leave the patch, so a macro-element assembles the exact Newton operator
with zero assembler surgery (no domain-level J̄ pass, no injected FE_Elements,
no ContactDomain-style machinery). Kernel reuse: `fbarScale2D` with the patch
J̄ and a generalized `addFbarCoupling2D` whose substitute gradient row is the
patch operator (patch-mean for T3, per-GP patch-P1 row for T6) instead of the
centroid g₀.

**Per-triangle formulation inside the macro:**
- **T3 pair (P4a):** patch-constant J̄ = v_patch/V_patch — dSNPO §15.1.9
  verbatim (their 2D recommendation IS 2-element patches; published anchors:
  Fig 15.8 localization, Fig 15.9 necking). Free-pair rank: exactly 3 RBM ✔.
- **T6 pair (P4b, by demand):** dilatation projected onto **P1 over the
  2-element patch** (L2, same quadrature, 6 GPs vs 3 P1 modes — a genuine
  least-squares smoothing). Free-pair rank: exactly 3 RBM, robust on
  distorted and conforming-curved pairs ✔. This lane is **novel/off-label**
  (equivalent mixed pair: P2 displacement / P1-per-macro pressure) → full
  adversarial gate + pressure battery when built.

**Findings that killed the alternatives (all numerically pinned):**

| # | Finding | Evidence |
|---|---|---|
| F1 | **dSNPO §15.1.9 is a disjoint-patch method, not nodal averaging** — the ADR's earlier "nodal" gloss conflated it with Bonet–Burton average-nodal-pressure. Disjointness is what makes the exact tangent patch-local (15.37/15.38). | book text pp. 665–669 |
| F2 | **F-bar-Patch (patch-constant J̄) does NOT cure the T6** — the P3 conformal modes re-center: u = a(z−z_p)², z_p = patch centroid, has ε_dev ≡ 0 pointwise and patch-mean dilatation 4Re(a(z̄_patch−z_p)) = 0 ⇒ zero energy under ANY constant-mean projection over ANY region. Free T6 pair: **5** zero modes (3 RBM + 2). Enlarging the averaging region can never fix the T6. | rank C4 = 5; energy injection E ≈ 1e-31 on distorted pairs |
| F3 | **Element-local P1 projection on the T6 is the identity** at the shipped 3-pt rule: #GP = dim(P1) = 3 and the interior points are P1-unisolvent, so the quadrature L2 projection interpolates ⇒ b̃ ≡ b at the GPs (max diff ~4e-15, straight AND curved — weights and curvature cancel algebraically). "‑formulation p1proj" would silently ship std. Corollary: **no element-local projection at this quadrature can yield any volumetric relief** — the cross-element coupling is the entire mechanism, not an implementation nuisance. | C6; twin reproduced incl. curved |
| F4 | **Patch-P1 over the pair restores full rank**: exactly 3 zero modes (symmetric, distorted, conforming-curved). Analytic proof (affine elements): ε_dev ≡ 0 piecewise + C0 on the shared edge ⇒ the two per-element analytic quadratics agree on a segment ⇒ identical (a degree-2 complex polynomial with >2 roots is null) ⇒ globally conformal ⇒ dilatation is a global linear field ⇒ reproduced exactly by patch-P1 ⇒ full volumetric energy. (Curved: result holds numerically; the written proof covers affine only.) | C5 = 3; twin distorted/curved = 3 |
| F5 | **Constraint counts** (T6 pair): std 6, patch-const 1 (over-relaxed → F2's 2 spurious modes), patch-P1 3 (= 1.5/element; mesh ratio r ≈ 2.67 vs std's locking 1.33). Heuristic only — the load-bearing evidence is F4 + F6. | C8 |
| F6 | **Inf-sup (T6 patch-P1)**: Dirichlet-pinned, mean-zero-pressure β_h **plateaus ≈ 0.46** for N = 2..5 structured pair-meshes, zero spurious pressure modes. Measured on straight structured meshes ONLY — no theoretical LBB proof exists for this pair; alternating-diagonal/unstructured/curved meshes unproven → P4b carries a pressure gate. | `p4_infsup_clean.py` |
| F7 | **patch-P2 is the same trap as F3 one level up**: 6 GPs = dim(P2) ⇒ identity ⇒ std ⇒ locks. Not a usable richer fallback. | twin ATTACK 2 |

**Rejected tangent strategies** (the original P4 shortlist, for the record):
1. *Nodal (overlapping) averaged dilatation* — cures both triangles but the
   exact tangent couples every vertex-sharing neighbor: needs domain-level
   injected FE_Elements (ADR-39-style) or an inexact tangent; κ-scaled dropped
   terms wreck Newton exactly where the cure matters; MP partitioning pain;
   plus the equal-order T3 pressure-oscillation record. Strictly dominated by
   the disjoint-patch macro, which gets exactness for free.
2. *Element-local P1disc projection (T6)* — dead three ways: F3 (identity at
   the shipped rule), no constraint relief even with a richer rule (still 3
   constraints/element), and P2/P1disc is the classically LBB-deficient pair.
3. *Lagged/staggered J̄* (freeze from last iteration) — destroys quadratic
   convergence and makes the consistent-tangent FD gate fail by design;
   against fork policy (no untestable mechanism ships).

**Implementation traps pinned for P4a/P4b** (from the adversarial gate):
- The K^(es) cross blocks are **stress-proportional** — they vanish at F = I,
  so the FD consistent-tangent gate MUST run at finite stress, and the element
  must declare an unsymmetric tangent (never symmetrize).
- **Conforming mesh assert**: a shared edge whose mid-nodes don't coincide
  silently cracks the patch and re-opens spurious modes (bit the twin's own
  curved test).
- **Patch partition policy**: leftover lone triangles = std (they lock);
  ≥3-element patches over-relax (F2 direction); user-facing form is a pair
  element (`element LadrunoXXX 4-or-9 nodes`) — apeGmsh emits natural pairs
  from split-quad meshing. Candidate ELE tag: **33021** (33019 pencilled H27,
  33020 SolidShell).
- F-bar power stays **½** (plane), as everywhere in this family.

Spike artifacts: `Ladruno_implementation/adr70_p4_spike/` — `p4_rank_spike.py`
(primary), `p4_adversarial.py` (independent twin: distorted pairs, energy
injection, 4-element over-relaxation trend, patch-P2 trap),
`p4_infsup_clean.py` (F6; NB its curved lane is degenerate — ignore, curved
inf-sup is genuinely open), `p4_curved_conforming.py` (F4 curved + the
cracked-mesh trap demo).
