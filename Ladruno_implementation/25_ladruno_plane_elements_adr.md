---
title: LadrunoCST & LadrunoQuad — unified 2D continuum elements
project: Ladruno
status: draft
priority: high
owner: nmora
tags:
  - implementation
  - element
  - solid
  - plane-stress
  - plane-strain
  - locking
  - 2d
---

# LadrunoCST & LadrunoQuad — unified 2D continuum elements

> Sibling document: [[26_ladruno_plane_frontier_adr]] — *why* this family exists
> beyond a port (the contribution thesis). **This** ADR is the Phase-0 substrate:
> the brick's proven element technology, lifted to 2D, byte-verified against
> upstream. Read 26 for the research direction; read 25 for what we build first.

## What

Two 2D continuum solid elements that mirror [[09_ladruno_brick]]'s design — the
**formulation is a parameter, not a class** — for plane-stress and plane-strain:

```
element('LadrunoQuad', tag, n1,n2,n3,n4, matTag, '-formulation', <std|bbar|ssp|eas>,
        '-geom', <linear|corot|finite>, '-thick', t, '-type', <PlaneStrain|PlaneStress>,
        '-hourglass', <stiffness|physical|viscous>)
element('LadrunoCST',  tag, n1,n2,n3,    matTag, '-geom', <linear|corot|finite>,
        '-thick', t, '-type', <PlaneStrain|PlaneStress>)
```

- **`LadrunoQuad`** — 4-node bilinear quad. `ELE_TAG_LadrunoQuad = 33007` (next
  in the element band after `LadrunoEmbeddedNode=33006`). Carries the full
  formulation × geometry menu, hourglass control, and the ASDConcrete crack-band
  `lch` handshake. This is where ~90% of the value lives.
- **`LadrunoCST`** — 3-node constant-strain triangle. `ELE_TAG_LadrunoCST = 33008`.
  Deliberately **thin**: `-formulation std` only (a 1-point triangle is
  rank-sufficient — no hourglass, and `bbar`/`ssp`/`eas` have nothing to average
  against), plus the geometry layer and material seam. It exists as the trivial
  baseline, the triangular-mesh fallback, and the future carrier for E-FEM
  ([[18_ladruno_efem_adr]]) — **not** as a contribution. See [[26_ladruno_plane_frontier_adr]] §CST
  for the honest verdict (plain CST is a dead end for fracture/plasticity).

This collapses the upstream 2D scatter — `Tri31`, `FourNodeQuad`, `SSPquad`,
`EnhancedQuad`, `bbarQuad`/`ConstantPressureVolumeQuad` — into two verifiable
classes, and adds the three things upstream has **none** of in 2D:
**corotational** plane solids, **finite-strain** plane solids, and the
**damage-scaled hourglass + crack-band** regularization.

**Scope — formulations (Quad):**

- `std`  — full 2×2 Gauss, displacement (reproduces upstream `FourNodeQuad`)
- `bbar` — mean-dilatation / selective volumetric integration (reproduces
  `ConstantPressureVolumeQuad` in plane strain; cures volumetric locking)
- `ssp`  — stabilized single-point (verbatim port of `UWelements/SSPquad::GetStab`;
  constant `Bnot`/`Kstab` condensed from the **initial** tangent — cures BOTH
  shear and volumetric locking across all ν; the explicit workhorse)
- `eas`  — Simo-Rifai enhanced assumed strain (4 internal modes Q1/E4, the 2D
  parent of the brick's E9; ported from `EnhancedQuad`). **Small-strain only in
  v1**; `eas`+`finite` is parser-reserved (Q1/E4 hourglasses in finite-strain
  compression — Wriggers & Reese 1996; see Risks).

**Scope — geometry (both elements), via [[solid_transformation_wrapper]] lifted to 2D:**

- `linear` — identity (byte-reduces to the direct small-strain kernel)
- `corot`  — element-level corotational; 2×2 polar / in-plane rotation
- `finite` — updated-Lagrangian; 2×2 **F**, Hencky log-strain via a new
  `LogStrain2D` material adaptor

**Scope — explicitly deferred:** everything in [[26_ladruno_plane_frontier_adr]]
(coupled-field gradient-damage / phase-field, VEM polygons, SBFEM, S-FEM). Those
are separate ADRs/PRs; this plan is the substrate they build on.

> **Design rule (inherited from the brick).** `-formulation` is a *single
> selector*, not independent booleans. B-bar ≡ selective integration of the
> volumetric term, so it cannot combine with `ssp` (also single-point) or be a
> flag on top of `std`. A selector makes invalid states unrepresentable.

## Why

1. **No corotational 2D solid exists in OpenSees.** `CorotCrdTransf` is
   frames-only; the 2D continuum elements (`Tri31`, `FourNodeQuad`, …) are all
   small-strain. Large-displacement plane analysis (membranes, snap-through
   panels, soil at finite rotation) has no home today.
2. **No finite-strain 2D solid exists.** *No* OpenSees plane element computes a
   2×2 deformation gradient **F**. Finite-strain plasticity / hyperelasticity in
   2D is simply unavailable — and 2D is where it should be developed and debugged
   before the brick (see point 5).
3. **No mesh-objective 2D concrete.** The brick's crack-band `lch` handshake +
   damage-scaled hourglass ([[11_brick_asdconcrete_integration.md]], PR #101)
   has no 2D counterpart. Plane RC walls / shear panels are the canonical 2D
   concrete problem.
4. **Scatter.** Choosing the right plane element today means knowing the folklore
   (`SSPquad` vs `EnhancedQuad` vs `ConstantPressureVolumeQuad` vs raw
   `FourNodeQuad`) and rebuilding the model to switch. One selector fixes that.
5. **2D is the proving ground for 3D.** This is the strategic reason. Meshes are
   cheap, a single Newton step is inspectable, and crack paths are *visible*.
   Every hard formulation the fork wants in the brick — finite-strain locking
   cures, and especially the coupled-field localization limiters of
   [[26_ladruno_plane_frontier_adr]] — is far cheaper to get right in 2D first,
   then lift. The plane family is the lab; the brick is production.

## Where

- **New code:**
  - `SRC/element/ladrunoPlane/LadrunoQuad.{cpp,h}` — the quad.
  - `SRC/element/ladrunoPlane/LadrunoCST.{cpp,h}` — the triangle.
  - `SRC/element/ladrunoPlane/OPS_LadrunoQuad.cpp`, `OPS_LadrunoCST.cpp` — parsers.
  - `SRC/element/ladrunoPlane/CMakeLists.txt`.
  - `SRC/material/nD/LogStrain2D.{cpp,h}` — 2D Hencky adaptor (2×2 **F**,
    plane-strain ε₃₃=0 / plane-stress σ₃₃=0). `ND_TAG_LogStrain2D = 33016`
    (RESERVED). May instead ship as `PlaneStrain`/`PlaneStress` getType views of
    the existing `LogStrainNDMaterial` (33010) — decided at P5 (see Risks).
- **Extend (the geometry layer):**
  - `SRC/element/solidTransformation/` — add 2D-aware corot/finite. The
    *interface* (`SolidTransformation.h`) already takes sized Vectors/Matrices
    with a runtime `numNodes`, so it is **not** 3D-locked. `SolidTransformationLinear`
    is pure identity → reuse as-is once the residual 3-hardcodes are dropped. The
    **corot** (3×3 polar) and **finite** (`invert3x3`, 3×3 **F**) *implementations*
    are 3D → add `SolidTransformation2DCorot` / `SolidTransformation2DFinite`
    (2×2 polar, 2×2 **F** + plane condensation). Factory `create(methodID, dim)`.
- **Reference (copy patterns from):**
  - `SRC/element/triangle/Tri31.{cpp,h}` — CST B-matrix, 1-GP, char-length.
  - `SRC/element/fourNodeQuad/FourNodeQuad.{cpp,h}` — 2×2 Gauss, std oracle.
  - `SRC/element/UWelements/SSPquad.{cpp,h}` — `GetStab` for `ssp`.
  - `SRC/element/fourNodeQuad/EnhancedQuad.{cpp,h}` — Q1/E4 EAS for `eas`.
  - `SRC/element/fourNodeQuad/ConstantPressureVolumeQuad.{cpp,h}` — `bbar` oracle.
  - `SRC/element/ladrunoBrick/LadrunoBrick.{cpp,h}` — the **architecture**:
    formulation dispatch, `isSinglePoint()`, damage-scaled `Kstab`, the three
    seams, sendSelf packing.
- **Ledgers / banner:** new rows in `LEDGER_implementations.md`
  (`LadrunoQuad`/`LadrunoCST`/`LogStrain2D`), `LEDGER_quirks.md` as we hit them,
  banner lines per shipped phase. Stamp the LADRUNO header on all new files
  (`Ladruno_scripts/stamp_headers.py`).

## How

### What ports vs. what is net-new

The brick mapping (verified by code review of `LadrunoBrick.cpp`):

| Reused (architecture / algorithm) | Net-new (2D numerics) |
|---|---|
| Formulation-enum dispatch switch | 3-component Voigt strain {11,22,12} (brick is 6) |
| `SolidTransformation` seams (sized) | B-matrix: 2 spatial gradients, 3 strain rows |
| SSP static-condensation logic | 1 hourglass mode {1,−1,1,−1} (brick has 4) |
| EAS inner-Newton condensation | Q1/E4 enhanced operator **M** (3×4) |
| `isSinglePoint()` material-eval economy | 2×2 (quad) / 1-pt (CST) quadrature |
| Tier-A damage-scaled hourglass | 2×2 **F** + plane condensation (finite) |
| sendSelf packing convention | char-length √A (quad) / √(2A) (CST) |
| std/bbar reduce-to-upstream gate pattern | `LogStrain2D` plane Hencky material |

**Voigt convention (state it once, use everywhere):** 2D ordering {11, 22, 12},
engineering shear γ₁₂ = 2ε₁₂. The strain–displacement operator per node *a*:

```
B_a = [ ∂N_a/∂x      0      ]      ε = B u  (Voigt, engineering shear)
      [    0      ∂N_a/∂y   ]
      [ ∂N_a/∂y   ∂N_a/∂x   ]
```

`bbar` replaces the volumetric rows with the element-mean gradient `∂N̄/∂x`,
`∂N̄/∂y`; shear row unchanged (the 2D analogue of the brick's `computeBbar`).

### Public API

As shown in **What**. Defaults: `-formulation std`, `-geom linear`,
`-type PlaneStrain`, `-thick 1.0`, no hourglass (only meaningful with a future
`uri`; for v1 `ssp` is the single-point path). `-type PlaneStress` requires a
material that supplies a plane-stress tangent (see material seam).

### Material seam (the 2D contract)

The element asks the nD material for a 2D view:

- **PlaneStrain:** `mat->getCopy("PlaneStrain")`. Works directly for materials
  with a plane-strain view (`LadrunoJ2` has it). A *3D-only* material
  (`ASDConcrete3D`, whose `getCopy("PlaneStress")` falls through to the base
  no-op — confirmed at `ASDConcrete3DMaterial.cpp:1849`) can still be driven in
  **plane strain** by feeding a 6-component strain with ε₃₃=γ₂₃=γ₁₃=0 and reading
  back the in-plane stress block — a thin "3D-in-plane-strain" adaptor inside the
  element. This is how plane-strain concrete works.
- **PlaneStress:** `mat->getCopy("PlaneStress")`. `LadrunoJ2` (nested ε₃₃ Newton
  → σ₃₃=0, already shipped) and `LadrunoRCConcrete` (PlaneStress view, 33015)
  supply it. Plain `ASDConcrete3D` does **not** → plane-stress concrete routes
  through `LadrunoRCConcrete` or a future plane-stress projection (Risks).
- **finite (`-geom finite`):** the element computes 2×2 **F**, hands it to
  `LogStrain2D` which forms the 2D Hencky strain, calls the inner plane material,
  and returns Cauchy σ + the spatial tangent. Plane-stress finite needs the σ₃₃=0
  thickness-stretch condensation (de Souza Neto Ch. 15 / Ch. 9).

### Hourglass & stabilization (Quad only)

A 4-node quad under 1-point integration has a **single** hourglass mode
γ = ¼(h − (h·x)bₓ − (h·y)b_y), h = {1,−1,1,−1} (Flanagan–Belytschko). The brick's
three hourglass flavours collapse in 2D:

- `ssp` = bbar + statically-condensed `Kstab` from the initial tangent (port
  `SSPquad::GetStab`). This is the v1 single-point element.
- **Damage-scaled `Kstab`** (the crack-band tie-in): degrade the constant elastic
  `Kstab` by s = max(1%, 1 − max(d_t, d_c)), damage read via a cached
  `setResponse("damage",…)` on the slot-0 material (null ⇒ s=1, zero change for
  elastic/J2). Verbatim the brick's Tier-A logic; `recvSelf` clears the cached
  Response so the raw material pointer can't dangle.

### Crack-band characteristic length

Override `getCharacteristicLength()`: **√A** for the quad (2×2-Gauss reference
area), **√(2A)** for the CST — the 2D analogues of the brick's ∛V and BezierTri6's
√(2A) ([[project_bezier_charlen]]). Replaces the `Element` base default
(min inter-node distance), which over-softens on graded/distorted meshes. Expose a
`charLength` response for verification.

### Geometry seam (2D SolidTransformation)

- **linear:** identity; the small-strain kernel routed through it must be
  bit-identical to the direct kernel (the seam proof).
- **corot:** `R` from the 2×2 polar decomposition of the in-plane deformation
  gradient (crib the closed-form 2×2 `F`→`atan2(f21−f12, f11+f22)` rotation
  extraction already in `ASDShellQ4CorotationalTransformation`); de-rotated small
  strain → reused small-strain plane material; exact residual + symmetric
  geometric tangent. No drilling DOF (continuum nodes have ndf=2).
- **finite:** updated-Lagrangian, 2×2 **F** = I + Σ uₐ ⊗ ∂Nₐ/∂X; `setTrialF` on
  `LogStrain2D`; assemble ∫ σ ∂N/∂x dv + the consistent spatial tangent. `bbar`+
  `finite` = 2D F-bar (F̄ = (J₀/J)^{1/2} F in 2D — note the **½** power, not ⅓;
  plane dilatation is 2D). Generally unsymmetric → `FullGeneral` advisory, as in
  the brick.

### Phasing (one PR each onto `ladruno`)

| Phase | Deliverable | Gate (oracle) |
|---|---|---|
| **P1** | `LadrunoQuad` std/bbar + `LadrunoCST` std, `-geom linear` | byte-reduce to `FourNodeQuad` / `ConstantPressureVolumeQuad` / `Tri31` ~1e-9; constant-strain patch; rank/3-rbm |
| **P2** | `ssp` + damage-scaled `Kstab` + crack-band `lch` (Quad) | `↔ SSPquad` ~1e-6 across ν∈{0,0.3,0.45,0.499}; mesh-objectivity (dissipation ∝ crack length); Tier-A scale == max(floor,1−d) |
| **P3** | `eas` (Quad, small strain) | distorted-mesh patch (∫M dV=0); reduce-to-std (α→0); bending-beats-std→Euler; near-incompressible non-locking |
| **P4** | `-geom corot` (both) | rigid-rotation >π objectivity; linear-bit-identical seam; FD-tangent-gap→0; load-driven cantilever vs finite <1% |
| **P5** | `-geom finite` + `LogStrain2D` (both) + 2D F-bar | FD consistent-tangent (unsym-aware); homogeneous-**F** patch vs oracle; reduce-to-linear at small strain; det F≤0 step-cut; volumetric-locking cantilever (std locks ν→0.5, F-bar doesn't) |

Test files: `tests/test_ladrunoQuad_{element,bbar,ssp,eas,corot,finite,asdconcrete}.py`,
`tests/test_ladrunoCST_{element,corot,finite}.py`, plus a Zone-B gmsh
`test_ladrunoQuad_asdconcrete_panel.py` (RC shear panel mesh-objectivity).

```cpp
// dispatch skeleton (mirrors LadrunoBrick)
enum class Formulation { STD, BBAR, SSP, EAS };
// getInitialStiff(): switch(formulation){ SSP→formSSP(); EAS→formEAS(); else std/bbar 2x2 }
// update(): if(isFinite()) updateFinite(); else if(SSP) ...; else std/bbar
// isSinglePoint(): SSP (and any future uri) → true → material eval'd once at slot 0
```

## Risks / open questions

> [!question]
> **Plane-stress concrete.** Plain `ASDConcrete3D` has no plane-stress view.
> P2/P5 plane-stress concrete must route through `LadrunoRCConcrete` (33015,
> PlaneStress) — or do we build a generic plane-stress *projection* wrapper
> (σ₃₃=0 static condensation of any 3D nD material, de Souza Neto Ch. 9)? The
> projection wrapper is more general and would also serve the brick's relatives.
> **Lean: ship plane-strain concrete in P2 (works via the 3D-in-plane-strain
> adaptor); defer plane-stress concrete to the projection-wrapper decision.**

> [!question]
> **`LogStrain2D` as a class vs. a view.** Make it a new nDMaterial (33016), or
> add `PlaneStrain`/`PlaneStress` getType views to `LogStrainNDMaterial` (33010,
> the `LadrunoJ2` dimensional-view pattern)? Views avoid a new tag and reuse the
> 3D return map; a separate class is cleaner for the 2×2-native **F** path.
> Decide at P5 after the geometry seam is wired.

> [!question]
> **EAS + finite hourglassing.** Simo-Rifai Q1/E4 develops spurious hourglass
> modes in finite-strain *compression* (Wriggers & Reese 1996; de Souza Neto
> Ch. 15) — exactly the concrete regime. v1 reserves `eas`+`finite`. The fix is
> stabilized-EAS or a Pian-Sumihara **assumed-stress hybrid** that does not
> hourglass — tracked as a frontier item in [[26_ladruno_plane_frontier_adr]] §T6.

- **2D F-bar power is ½, not ⅓.** Plane dilatation J = det(2×2 F); the F-bar
  scaling is (J₀/J)^{1/2}. Copy-pasting the brick's ^{1/3} is a bug. Pin it with
  a near-incompressible plane-strain cantilever.
- **CST is honestly low-value.** It will *pass* its gates (it reduces to `Tri31`)
  but it locks volumetrically and biases localization to mesh lines. Ship it as
  the baseline/E-FEM-carrier; do not over-invest. See [[26_ladruno_plane_frontier_adr]].
- **ndf consistency.** Continuum nodes are ndf=2; corot adds **no** drilling DOF
  (unlike the shell transformations, which assume ndf=6 — do **not** reuse
  `ASDEICR`'s 6-DOF machinery, only its 2×2 `F` rotation extraction).
- **Backwards compatibility:** new elements, new tags — no existing-model impact.

## Implementation log

*(filled in as phases land; move to `Ladruno_internal/` when complete)*
