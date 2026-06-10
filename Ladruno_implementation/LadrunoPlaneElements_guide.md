---
title: "LadrunoQuad & LadrunoCST — 2D Plane Continuum Elements"
project: Ladruno
type: reference / user guide
status: v1 shipped (ADR 25 Phase 1 std/bbar + Phase 2 ssp + crack-band lch + Tier-A damage; 2D linear kinematics). LadrunoQuad 33007, LadrunoCST 33008
element: "element LadrunoQuad / element LadrunoCST"
classTag: 33007 (LadrunoQuad), 33008 (LadrunoCST)
related:
  - "[[25_ladruno_plane_elements_adr]]"
  - "[[26_ladruno_plane_frontier_adr]]"
  - "[[ladruno_continuum_elements_guide]]"
  - "[[09_ladruno_brick]]"
  - "[[solid_transformation_wrapper]]"
  - "[[LEDGER_implementations]]"
tags:
  - guide
  - element
  - plane-stress
  - plane-strain
  - quad
  - cst
  - triangle
  - bbar
  - ssp
  - crack-band
  - 2d-continuum
updated: 2026-06-10
---

# LadrunoQuad & LadrunoCST — 2D plane continuum elements

`element LadrunoQuad` (classTag **33007**) is a unified **4-node bilinear** plane
(plane-stress / plane-strain) continuum element — the 2D sibling of
[[09_ladruno_brick|LadrunoBrick]], where the **formulation is a parameter, not a
class**: one element name spans standard displacement, B-bar, and stabilized
single-point. `element LadrunoCST` (classTag **33008**) is its **thin 3-node
constant-strain-triangle** sibling: `std` only (a 1-point triangle is
rank-sufficient — there is nothing to stabilize), shipped as the baseline /
triangular-mesh fallback / future E-FEM carrier.

Both are **small-strain, geometrically linear** in v1 (ADR 25 Phase 1/2). They are
careful ports of the upstream plane elements and reduce to them bit-for-bit:

| Ladruno | `-formulation` | reduces to (upstream) | to |
|---|---|---|---|
| `LadrunoQuad` | `std` | `FourNodeQuad` (cmd `quad`) | ~1e-9 |
| `LadrunoQuad` | `bbar` | `ConstantPressureVolumeQuad` | (mean-dilatation) |
| `LadrunoQuad` | `ssp` | `SSPquad` (UWelements) | ~1e-6 |
| `LadrunoCST` | (std) | `Tri31` | ~1e-9 |

This is the single user/developer reference: **theory** (kinematics, the three
formulations, the crack-band length), **capabilities** (the option grammar),
**modeling guidance** (which formulation, plane-stress vs strain, mass, softening,
initial strain), **implementation**, **validation**, and **use cases**. The design
record is [[25_ladruno_plane_elements_adr|ADR 25]] (basic port) with the
research/contribution frontier in [[26_ladruno_plane_frontier_adr|ADR 26]]. For
choosing *between* continuum elements see [[ladruno_continuum_elements_guide]].

> [!important] CST is honestly low-value — reach for it deliberately
> A plain constant-strain triangle **volumetric-locks** (it cannot represent the
> isochoric modes a near-incompressible material needs) and **mesh-biases
> localization**. Use `LadrunoCST` for triangular-mesh regions you can't avoid, as
> a coarse baseline, or as the E-FEM carrier — not as a workhorse. For real 2D
> work prefer `LadrunoQuad` (quads) or `BezierTri6` (quadratic triangle). See
> [[26_ladruno_plane_frontier_adr]] §CST.

---

## 1. Theory

### 1.1 Kinematics — the shared small-strain core

Both elements work in a `ndm = 2, ndf = 2` model (the parser **refuses** anything
else). DOFs are the in-plane translations `(u_x, u_y)` at each node — 8 for the
quad, 6 for the CST. The strain is the engineering Voigt triple
`ε = [ε_xx, ε_yy, γ_xy]ᵀ = B u`, the stress `σ = [σ_xx, σ_yy, τ_xy]ᵀ` comes from the
assigned `nDMaterial`, and the element assembles

```
residual   P = ∫ Bᵀ σ dV − f_body − f_pressure − Q       (getResistingForce)
tangent    K = ∫ Bᵀ D B dV                                (getTangentStiff,  D = material tangent)
```

over the integration rule. **Quad** uses the standard `2×2` Gauss rule (4 points)
for `std`/`bbar` and a **single centroid** evaluation for `ssp`. **CST** uses one
centroid point (the strain is constant over the element — the defining CST
property). `dV = detJ · w · t`, with `t` the out-of-plane **thickness** (`-thick`,
default 1.0).

### 1.2 The three quad formulations

**`std` — standard displacement (default).** Full `2×2` Gauss, bilinear shape
functions, `B` straight from the shape-function gradients. Reduces to upstream
`FourNodeQuad`. The honest default; **volumetric-locks** as `ν → 0.5` in plane
strain and bends stiffly on coarse meshes.

**`bbar` — mean-dilatation B-bar (plane strain only).** The dilatational part of
`B` is replaced by its element average, curing volumetric locking while preserving
constant-strain completeness (the patch test still passes). The 2D dilatation
factor is **½** (not the brick's ⅓): for node `a`,

```
c_x = ½(b̄_x − b_x),   c_y = ½(b̄_y − b_y)        (b̄ = element-mean gradient)
B(ε_xx) row:  [ b_x + c_x ,  c_y ]
B(ε_yy) row:  [ c_x ,  b_y + c_y ]
B(γ_xy) row:  [ b_y ,  b_x ]                       (shear unchanged)
```

`bbar` is **refused with `-type PlaneStress`** at the parser — volumetric locking is
a plane-strain / incompressible phenomenon, so B-bar there is meaningless.

**`ssp` — stabilized single-point.** A verbatim port of `UWelements/SSPquad::GetStab`:
the material is evaluated **once** at the centroid (`K = 4 J₀ t · Mᵀ C M`), and a
statically-condensed **hourglass stabilization** `K_stab` — built from the
**initial** material tangent, so it is geometry-fixed and computed once in
`setDomain` — suppresses the spurious mode. `ssp` cures **both** shear and
volumetric locking across all `ν` and is the cheapest quad for nonlinear /
explicit work (one material point instead of four). It is the 2D analogue of
`LadrunoBrick -formulation ssp`.

> [!note] `eas` is reserved
> `-formulation eas` (Simo–Rifai enhanced assumed strain, the 2D Q1/E4 analogue of
> the brick's EAS) is **parsed but refused** with a "not yet implemented" hint
> (ADR 25 Phase 3). The API slot is reserved so it drops in without a grammar change.

### 1.3 CST — the constant-strain triangle

`LadrunoCST` has linear (area-coordinate) shape functions, so `B` is **constant**
over the element and the single centroid point is exact. There is **no formulation
axis** — `bbar`/`ssp`/`eas` would have nothing to average or stabilize on a
rank-sufficient 1-point triangle. It reduces to upstream `Tri31`.

### 1.4 Crack-band characteristic length

Both override `getCharacteristicLength()` so a regularized softening material
(ASDConcrete crack-band, Lemaitre, …) gets a **geometry-true** length, not the
base class's edge-distance default:

```
LadrunoQuad:  lch = √A          (A = ∫ detJ dξ — the true element area)
LadrunoCST:   lch = √(2A)       (matches the BezierTri6 triangle convention)
```

Exposed as the `charLength` response (§7) and consumed automatically by crack-band
materials for mesh-objective dissipation.

### 1.5 Tier-A damage-scaled stabilization (ssp + softening)

Because the `ssp` `K_stab` is built from the **elastic** initial tangent, a fully
cracked element would keep spurious elastic hourglass control. Tier-A degrades it:

```
K_stab,eff = s · K_stab,elastic,   s = max(0.01, 1 − max(d_t, d_c))
```

`d_t/d_c` are read from a cached `"damage"` response probe on the centroid material
(materials with no damage channel return nothing → `s = 1`, so elastic/J2 are
**unchanged**). The 1% floor keeps a sliver of hourglass control in a fully-damaged
element. This mirrors [[09_ladruno_brick|LadrunoBrick]]'s Tier-A and fires only for
the `ssp` formulation.

---

## 2. Capabilities — the option grammar

```
element LadrunoQuad $tag $n1 $n2 $n3 $n4 $matTag
        [-formulation {std | bbar | ssp | eas}]      # default std; eas reserved
        [-type {PlaneStrain | PlaneStress}]          # default PlaneStrain
        [-thick $t] [-rho $r] [-body $bx $by] [-pressure $p]

element LadrunoCST  $tag $n1 $n2 $n3 $matTag
        [-type {PlaneStrain | PlaneStress}]
        [-thick $t] [-rho $r] [-body $bx $by] [-pressure $p]
```

| Token | Meaning | Default |
|---|---|---|
| `$n1..n4` / `$n1..n3` | corner nodes, **CCW**, each `ndf = 2` | — |
| `$matTag` | an `nDMaterial` that supports the requested plane view (`getCopy("PlaneStrain"/"PlaneStress")`) | — |
| `-formulation` | quad only: `std`/`bbar`/`ssp` (`eas` reserved) | `std` |
| `-type` | `PlaneStrain` or `PlaneStress` | `PlaneStrain` |
| `-thick $t` | out-of-plane thickness | `1.0` |
| `-rho $r` | mass density override (else the material's `getRho()`) | `0.0` |
| `-body $bx $by` | constant body-force per unit volume | `0, 0` |
| `-pressure $p` | uniform edge pressure (consistent nodal load) | `0.0` |

The model **must** be `ndm 2, ndf 2`; nodes must carry exactly 2 DOFs or
`setDomain` refuses. `bbar` + `PlaneStress` is refused; `eas` is refused. Unknown
options are **warn-and-ignore** (permissive OpenSees style) — watch for a typo'd
flag silently building a default element. Build nodes CCW so `detJ > 0` (see §6).

---

## 3. Choosing the formulation (quad)

| Situation | Use | Why |
|---|---|---|
| General-purpose, well-shaped mesh | **`std`** | the honest baseline; reduces to FourNodeQuad |
| Near-incompressible plane strain (`ν → 0.5`, J2 plasticity, undrained soil) | **`bbar`** or **`ssp`** | cure volumetric locking; `std` locks |
| Coarse bending-dominated mesh | **`ssp`** | stabilized; far better coarse-bending than `std` |
| Nonlinear / explicit (expensive material) | **`ssp`** | **one** material eval/element instead of four |
| Softening / crack-band concrete | **`ssp`** (+ a crack-band material) | single-point + Tier-A damage-scaled `lch` handshake |
| You need a triangle | `LadrunoCST` or `BezierTri6` | CST locks/biases — prefer `BezierTri6` if you can |

`bbar` and `ssp` both relieve volumetric locking; `ssp` additionally relieves shear
(coarse-bending) locking and is cheaper, so for nonlinear work `ssp` is usually the
better pick. `bbar` stays closest to the textbook B-bar element for verification.

---

## 4. Plane stress vs plane strain

`-type` selects the material view via `material->getCopy(type)`, so the **material**
enforces the constraint (`σ_zz = 0` for plane stress; `ε_zz = 0` for plane strain) —
the element is view-agnostic. Make sure the assigned `nDMaterial` actually offers
the requested view (`ElasticIsotropic`, J2, ASDConcrete plane views, …); a material
that returns 0 for `getCopy("PlaneStrain")` will fail at construction.

- **Plane strain** — long prismatic bodies (dams, tunnels, retaining walls, a slice
  of a long structure). The locking regime; reach for `bbar`/`ssp`.
- **Plane stress** — thin in-plane-loaded plates / shear walls / membranes. No
  volumetric locking, so `bbar` is disallowed; use `std` or `ssp`.

---

## 5. Cross-cutting modeling

### 5.1 Mass & explicit dynamics
`getMass()` builds a **lumped diagonal** mass (shape-function-weighted), **not** the
consistent mass upstream `FourNodeQuad` uses. Total mass is conserved and the
diagonal form is `DiagonalSOE`-safe for explicit central-difference. Density comes
from `-rho` or the material's `getRho()`; zero density ⇒ zero mass matrix.
**Consequence:** "reduces to FourNodeQuad" holds for **statics** (stiffness /
stress / force) — the dynamic mass differs by design. For `ssp` the lumped mass
uses the SSP `J₀ + J₁ξ + J₂η` weights.

### 5.2 Near-incompressibility
Use `bbar` (plane strain) or `ssp`. `std` plane-strain displacements stiffen
sharply as `ν → 0.5`; the validation suite pins `bbar` strictly more flexible than
`std` at `ν = 0.4999`.

### 5.3 Softening / concrete
Use `ssp` with a crack-band material. The `lch = √A` handshake (§1.4) makes
dissipation mesh-objective; Tier-A (§1.5) keeps hourglass control as the element
cracks. CST softening is **not recommended** (mesh-biased localization).

### 5.4 Initial strain / prestrain — at the *material* level
**Neither element has an element-level initial-strain option** (no `-initStrain`;
they feed the full `B u` straight to `setTrialStrain`). Apply prestrain / eigenstrain
by **wrapping the `nDMaterial`**:

- **`InitStrainNDMaterial`** — a constant initial strain `ε₀` subtracted inside the
  material. Now **dimension-general** (`PlaneStrain` / `AxiSymmetric` views were
  added in the fork, [[LEDGER_implementations]] / PR #222), so it composes with
  these plane elements in 2D, not just 3D.
- **`StagedStrain`** (classTag 33014, [[project_initdefgrad_staged]]) — stress-free
  staged activation (`ε_rel = ε − ε₀`, auto-captured at activation), and it composes
  with `InitStrain` for a prestrained-and-staged member.

```python
ops.nDMaterial('ElasticIsotropic', 1, E, nu)
ops.nDMaterial('InitStrain', 2, 1, eps0_xx, eps0_yy, eps0_xy)   # wrap with ε₀ (2D-general)
ops.element('LadrunoQuad', 1, 1,2,3,4, 2, '-type', 'PlaneStrain')
```

### 5.5 Recording
Per-element `force`, per-Gauss-point `stresses`/`strains`, `charLength`, and the
per-point `material`/`integrPoint` response (§7).

---

## 6. Gotchas & limits

- **Node ordering must be CCW** (`detJ > 0`). `shapeFunction()` does **not** guard
  `detJ ≤ 0` (matching upstream `Tri31`/`FourNodeQuad`): a clockwise or degenerate
  element silently produces a wrong (possibly indefinite) stiffness. `getCharacteristicLength`
  *does* guard a non-positive area and falls back. If a model won't converge on a
  triangular/quad patch, check winding first.
- **`bbar` is plane-strain only**; **`eas` is reserved** (both refused at parse).
- **CST has no formulation axis** — passing `-formulation` to `LadrunoCST` is a
  warn-and-ignored unknown option (the parser has no such flag).
- **Mass is lumped, not consistent** (§5.1) — a dynamic cross-check against
  `FourNodeQuad` will differ.
- **Small strain only** in v1 — no `-geom corot/finite` yet (§8).
- **Material must support the requested plane view** or construction fails.

---

## 7. Diagnostics & responses

`eleResponse $tag <name>` / `recorder Element -ele $tag -<name>`:

| Response | Aliases | Size (Quad / CST) | Meaning |
|---|---|---|---|
| `force` | `forces` | 8 / 6 | element resisting force vector |
| `stresses` | `stress` | 12 / 3 | per-GP stress `[σxx,σyy,τxy]` (Quad: 4 GP; `ssp` mirrors the centroid to all 4) |
| `strains` | `strain` | 12 / 3 | per-GP strain (same layout) |
| `charLength` | `characteristicLength` | 1 | crack-band length (`√A` quad, `√(2A)` CST) |
| `material $i …` | `integrPoint $i …` | — | forward to integration-point `i`'s material response (`ssp` always uses point 1) |

`setParameter "pressure" $p` (live edge-pressure update) and per-point
`material $i <param>` are supported for parametric studies.

---

## 8. Deferred / not-yet-built

| Feature | Status |
|---|---|
| `-geom corot` / `-geom finite` (large displacement / finite strain) | deferred — drops in via the 2D [[solid_transformation_wrapper]] seam (ADR 25 P4/P5: `SolidTransformation2DCorot` / `2DFinite`, F-bar power **½**) |
| `-formulation eas` (Simo–Rifai Q1/E4) | reserved, parser-refused (ADR 25 P3) |
| `LogStrain2D` finite-strain material adaptor (`33016`) | draft / reserved (may ship as `PlaneStrain`/`PlaneStress` views of `LogStrainNDMaterial 33010`) |
| Gradient-damage / phase-field localization (T1/T2) | research frontier — [[26_ladruno_plane_frontier_adr]] |
| Element-level `-initStrain` | not planned — use `InitStrainNDMaterial` / `StagedStrain` (§5.4) |

---

## 9. Implementation map

| Concern | Where |
|---|---|
| Elements | `SRC/element/ladrunoPlane/LadrunoQuad.{cpp,h}`, `LadrunoCST.{cpp,h}` |
| Parsers | `OPS_LadrunoQuad.cpp`, `OPS_LadrunoCST.cpp` (flag-based; reject `eas`, `bbar`+PlaneStress, non-2D models) |
| Shape fns / `B` | `shapeFunction()` (bilinear quad / linear CST), `formB()` (quad std/bbar), `computeShapeBar()` (quad B-bar mean gradients) |
| SSP | `computeSSP()` (port of `SSPquad::GetStab`: `J₀/J₁/J₂`, centroid `Mmem`, condensed `Kstab` from the **initial** tangent) — built once in `setDomain` |
| Crack-band | `getCharacteristicLength()` (`√A` / `√(2A)`); `damageScale()` reads the cached `"damage"` probe (ssp) |
| Residual / tangent | `getResistingForce` (`∫Bᵀσ` − body − pressure − Q), `getTangentStiff` (`∫BᵀDB`), `getInitialStiff` (cached `*Ki`) |
| Mass / loads | `getMass` (lumped), `addLoad` (`SelfWeight`), `setPressureLoadAtNodes` (consistent edge pressure) |
| Serialization | `sendSelf`/`recvSelf` carry tag/thickness/body/pressure/rho/formulation/planeType + Rayleigh + material tags; SSP state & `B` recomputed in `setDomain` on receive (the cached damage probe is dropped in `recvSelf` so it can't dangle) |
| Plumbing | `classTags.h` (33007 / 33008), `FEM_ObjectBrokerAllClasses.cpp`, `OpenSeesElementCommands.cpp`, CMake — see [[LEDGER_vanilla_files]] |

---

## 10. Validation

Zone-A battery `tests/test_ladrunoPlane_element.py`:

- **Headline reductions:** `std` ↔ `FourNodeQuad` (plane strain + plane stress,
  ~1e-9 disp/stress/force); `LadrunoCST` ↔ `Tri31`; default formulation == `std`.
- **`ssp` ↔ `SSPquad`** to ~1e-6 across `ν ∈ {0, 0.3, 0.45, 0.499}` plane strain +
  plane stress (the stabilization is the validated quantity — patch/rank can't
  prove it, the cross-check can).
- **Rank / zero-energy:** quad (`std`, `bbar`, `ssp`) and CST each have exactly 3
  rigid-body modes — no spurious mechanisms.
- **Volumetric-locking relief:** at `ν = 0.4999`, `bbar` is strictly more flexible
  than `std`.
- **Constant-stress patch** (`std` + `bbar`): every Gauss point reports the
  closed-form uniaxial stress (B-bar preserves completeness).
- **Crack-band `lch`:** `√(a·b)` for a stretched quad, `√(p·q)` for a right-triangle
  CST (≠ min edge).
- **Guards:** `eas` refused, `bbar`+PlaneStress refused, `ssp` builds.

---

## 11. Use cases & recipes

```python
# ── plane-strain, near-incompressible, implicit (dam slice / undrained clay) ──
ops.model('basic', '-ndm', 2, '-ndf', 2)
# ... CCW nodes ...
ops.nDMaterial('ElasticIsotropic', 1, E, 0.49)
ops.element('LadrunoQuad', 1, 1,2,3,4, 1, '-type', 'PlaneStrain', '-formulation', 'bbar')

# ── thin shear wall, plane stress ──
ops.element('LadrunoQuad', 2, 5,6,7,8, 1, '-type', 'PlaneStress', '-thick', 0.2)

# ── nonlinear / softening concrete panel: ssp + crack-band material ──
ops.nDMaterial('ASDConcrete3D', 10, ...)            # a plane view that supports "damage"
ops.element('LadrunoQuad', 3, 9,10,11,12, 10, '-type', 'PlaneStress', '-formulation', 'ssp')
#   → single material eval, √A crack-band lch, Tier-A damage-scaled hourglass

# ── triangular-mesh region (use sparingly) ──
ops.element('LadrunoCST', 4, 13,14,15, 1, '-type', 'PlaneStrain', '-thick', 0.5)

# ── prestrain via the material wrapper (no element flag) ──
ops.nDMaterial('InitStrain', 20, 1, 1.0e-3, 1.0e-3, 0.0)
ops.element('LadrunoQuad', 5, 16,17,18,19, 20, '-type', 'PlaneStrain')

# ── self-weight / body force ──
ops.element('LadrunoQuad', 6, 20,21,22,23, 1, '-body', 0.0, -9.81*rho)
```

### apeGmsh
For a 2D structural mesh, emit `LadrunoQuad` per quad face (`-formulation ssp` for
nonlinear/softening, `bbar` for verification plane-strain), `LadrunoCST` only where
the mesh is genuinely triangular, passing `-type`, `-thick`, and the resolved
`nDMaterial` tag.

---

## 12. References & related
- **[[25_ladruno_plane_elements_adr|ADR 25]]** — the basic-element port (formulations,
  the 2D B-bar/SSP derivations, the geometry-seam plan).
- **[[26_ladruno_plane_frontier_adr|ADR 26]]** — the research frontier (gradient-damage,
  phase-field, VEM, SBFEM, S-FEM) and the honest CST assessment.
- **[[09_ladruno_brick|LadrunoBrick]]** — the 3D sibling these mirror (same
  `-formulation` philosophy, SSP, Tier-A, crack-band `lch`).
- **[[ladruno_continuum_elements_guide]]** — choosing *between* continuum elements
  (Quad / CST / BezierTri6 / BezierTet10 / Brick).
- **[[solid_transformation_wrapper]]** — the corot/finite geometry seam the plane
  elements will consume (ADR 25 P4/P5).
- **Upstream:** `FourNodeQuad` · `ConstantPressureVolumeQuad` · `SSPquad` · `Tri31`.
</content>
