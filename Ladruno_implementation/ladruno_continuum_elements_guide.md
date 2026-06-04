---
title: Ladruno Continuum Elements — Modeling & FE-Selection Guide
project: Ladruno
status: living-reference
priority: high
owner: nmora
tags:
  - reference
  - guide
  - element
  - selection
  - solid
  - continuum
  - locking
  - explicit
  - finite-strain
aliases:
  - "Continuum element guide"
  - "Element selection guide"
  - "FE selection"
---

# Ladruno Continuum Elements — Modeling & FE-Selection Guide

> [!abstract] What this document is
> A **single decision desk** for the fork's three displacement continuum
> elements — `BezierTri6` (2D), `BezierTet10` (3D tet) and `LadrunoBrick`
> (3D hex). It answers the modeling questions *before* you reach a per-element
> reference: **which element, which formulation, which geometry method, which
> mass, which solver**, for the analysis in front of you. The deep theory for
> each element lives in its own doc — this guide stays at the altitude of
> *intended use* and *selection procedure* and links down.
>
> **Per-element deep docs (the references this guide points to):**
> [[LadrunoBrick_reference]] · [[04_bezier_elements|BezierTri6 ADR]] ·
> [[06_bezier_tet10|BezierTet10 ADR]] · [[14_bezier_tet10_corot|BezierTet10 corot]] ·
> [[11_brick_asdconcrete_integration|Brick + concrete/softening]] ·
> [[solid_transformation_wrapper|geometry seams]] ·
> [[09_finite_strain_material_wrapper|finite-strain materials]].

---

## Table of contents

1. [[#1 · The roster at a glance]]
2. [[#2 · The selection axes]]
3. [[#3 · The decision procedure]]
4. [[#4 · Intended-use profiles (per element)]]
5. [[#5 · Cross-cutting modeling guidance]]
6. [[#6 · Quick-start recipes]]
7. [[#7 · Anti-patterns & common mistakes]]
8. [[#8 · References & maintenance]]

---

## 1 · The roster at a glance

All three are **pure-displacement continuum elements** that delegate the
constitutive law to an `NDMaterial`. They share the fork's design DNA:
non-negative lumped mass for explicit dynamics, geometry-true characteristic
length for crack-band regularization, and the Ladruno recorder contract.

| | **BezierTri6** | **BezierTet10** | **LadrunoBrick** |
|---|---|---|---|
| **classTag** | `33000` | `33001` | `33002` |
| **Dim / nodes / DOF** | 2D / 6 / 12 | 3D / 10 / 30 | 3D / 8 / 24 |
| **Topology** | quadratic triangle | quadratic tetrahedron | trilinear hexahedron |
| **Order** | quadratic (Bernstein) | quadratic (Bernstein) | linear |
| **Material type** | `PlaneStress`/`PlaneStrain` | `ThreeDimensional` | `ThreeDimensional` |
| **Anti-locking** | `-bbar` (plane strain) | `-bbar` (+ `-fbar`) | `-formulation bbar`/`uri`/`ssp`/`eas` |
| **Geometry methods** | linear only | `linear`/`corot`/`finite` | `linear`/`corot`/`finite` |
| **Default mass** | lumped `ρVe/6` (all +) | lumped `ρVe/10` (all +) | consistent (`-lumped` for diag) |
| **Explicit-ready** | ✅ (lumped is +) | ✅ (lumped is +) | ✅ (`uri`/`ssp` + `-lumped`) |
| **Quadrature** | fixed 3-pt | fixed 4-pt | 2×2×2 (or 1-pt for `uri`/`ssp`) |
| **Char. length** | `sqrt(2·A)` | `cbrt(6·V)` | `cbrt(V)` |
| **Curved edges** | ❌ straight-sided v1 | ❌ straight-sided v1 | n/a |
| **Can host embedded rebar** | ✅ | ✅ | ✅ (see [[20_ladruno_embedded_reinforcement_adr]]) |
| **Deep doc** | [[04_bezier_elements]] | [[06_bezier_tet10]] | [[LadrunoBrick_reference]] |

> [!note] One-line identity
> - **`BezierTri6`** — the quadratic 2D workhorse: smooth stress fields and
>   explicit-ready lumped mass from standard T6 meshes.
> - **`BezierTet10`** — the same idea in 3D, for geometry that is painful to
>   hex-mesh; now with `corot`/`finite` kinematics + F-bar.
> - **`LadrunoBrick`** — the *unified* hex: one element, four formulations
>   (`std`/`bbar`/`uri`/`ssp`/`eas`) × three geometry methods, bit-identical to
>   the upstream `Brick`/`bbarBrick`/`SSPbrick` where they overlap.

---

## 2 · The selection axes

Element choice is **not one decision** — it factors into a handful of nearly
independent axes. Resolve each, and the element + flags fall out. The axes, in
roughly the order they bind:

### Axis A — Dimensionality
- **2D (plane stress / plane strain / membrane slice)** → **`BezierTri6`**.
  It is the only 2D continuum element in the fork's set.
- **3D** → `BezierTet10` **or** `LadrunoBrick` — decided by Axis B.

### Axis B — Mesh type / geometry (the 3D fork)
This is the dominant 3D decision and it is driven by **meshability**, not by
the physics:
- **Blocky / sweepable / structured geometry** (prisms, walls, dams, regular
  solids) → **`LadrunoBrick`** (hex). A hex mesh gives more accuracy per DOF and
  per-element cost than tets, and the brick's formulation menu is richer.
- **Complex / organic / CAD-import geometry** that does not hex-mesh cleanly
  (re-entrant corners, fillets, arbitrary STEP solids) → **`BezierTet10`**.
  Tet meshing is automatic; the quadratic Bézier basis buys back much of the
  accuracy that linear tets (`FourNodeTetrahedron`) throw away.
- **Mixed** → mesh the bulk with hex and transition regions with tets; both
  carry the same recorder contract and char-length convention.

> [!tip] Hex if you can, tet if you must
> A quadratic tet (Tet10) is far better than a linear tet, but a well-shaped hex
> still wins on accuracy-per-DOF for the same problem. Reach for `BezierTet10`
> when the **geometry** forces your hand, not as a default 3D element.

### Axis C — Element order / accuracy
- **Linear hex** (`LadrunoBrick`) — 1 quadratic-accurate displacement field per
  3 edges; needs the anti-locking formulations to behave in bending and near
  incompressibility.
- **Quadratic tet/tri** (`Bezier*`) — the displacement field is quadratic, so
  bending and stress recovery are **naturally good** without enhanced strains;
  the price is more nodes/DOF per element and (today) straight-sided geometry.

### Axis D — Analysis type (implicit vs explicit)
- **Implicit** (static push-over, quasi-static, implicit dynamics) — any element;
  prefer well-conditioned tangents. For the brick: `bbar`/`std`/`ssp`/`eas`.
- **Explicit** (impact, blast, contact-rich, snap-back fallback) — needs a
  **non-negative lumped mass** and ideally a cheap internal force:
  - `BezierTri6` / `BezierTet10` are *designed* for this — their Bernstein
    lumped mass (`ρVe/6`, `ρVe/10`) is strictly positive (the headline reason
    these elements exist). Add `-lumped` is **not** needed; lumped is the
    default — pass `-cMass` only to force consistent mass.
  - `LadrunoBrick` → `-formulation uri -hourglass physical` (8 damage points) or
    `-formulation ssp` (cheapest, 1 material eval), **plus `-lumped`**.
  - Pair with an explicit integrator ([[central_difference_ladruno_guide|CentralDifferenceLadruno]]
    / [[explicit_bathe_integrators_guide|ExplicitBathe]]) and a `criticalTimeStep()`-driven `dt`.

### Axis E — Locking regime
- **Near-incompressible** (`ν → 0.5`, undrained soil, J2 at the isochoric yield
  surface, rubber) → use a B-bar / F-bar / SSP cure on **every** element:
  - `BezierTri6 -bbar` (plane strain only — plane stress has no volumetric
    constraint, so `-bbar` is rejected there).
  - `BezierTet10 -bbar`.
  - `LadrunoBrick -formulation bbar` (or `ssp` across all ν, well-conditioned).
- **Bending-dominated, coarse mesh** → quadratic Bézier elements handle it
  naturally; for the linear brick use `ssp` or `eas` (true Simo–Rifai enhanced
  assumed strain), **not** `bbar` (which cures volumetric locking only).
- **Both** (thin near-incompressible bending) → `LadrunoBrick -formulation ssp`,
  or refine a quadratic Bézier mesh.

### Axis F — Kinematic regime
- **Small strain, small rotation** → `-geom linear` (default everywhere).
- **Large rotation, small strain** (rotating member, buckling of a stiff part)
  → `-geom corot` — available on `BezierTet10` and `LadrunoBrick` (**not**
  `BezierTri6`, which is small-displacement only). Small-strain materials are OK
  as long as strains stay small.
- **Large strain** (forming, soft soil, hyperelastic, finite plasticity)
  → `-geom finite` — `BezierTet10` and `LadrunoBrick`. Requires a
  **`FiniteStrainNDMaterial`** (e.g. `nDMaterial LogStrain` wrapping a 3D
  material; see [[09_finite_strain_material_wrapper]]). `-geom finite -bbar`
  selects the **F-bar** element for near-incompressible large strain — its
  tangent is **generally unsymmetric**, so use an unsymmetric solver
  (`system FullGeneral`/`UmfPack`/`SparseGEN`).

### Axis G — Material / softening
- **Linear-elastic / hardening plasticity** → any element, any formulation.
- **Softening (concrete, ductile damage)** → crack-band regularization needs a
  geometry-true characteristic length; all three elements override
  `getCharacteristicLength()` so the material's `lch` is mesh-true. For the
  brick + `ASDConcrete3D` the single-point `Kstab` **damage-scales** (Tier-A) so
  `ssp`/`uri` are usable under softening — but full integration (`std`/`bbar`)
  stays the most robust. **Monitor `hourglassEnergy` < 5–10 %.** See
  [[11_brick_asdconcrete_integration]].

---

## 3 · The decision procedure

Walk the axes top-to-bottom. The first three pick the **element**; the rest pick
the **flags**.

```text
START
 │
 ├─ 2D problem? ───────────────────────────► BezierTri6
 │        └─ flags: PlaneStress | PlaneStrain ; -bbar if PlaneStrain & ν→0.5
 │                  default lumped mass (explicit-ready); -cMass for implicit accuracy
 │
 └─ 3D problem
      │
      ├─ Geometry hex-meshable (blocky/sweepable)? ──► LadrunoBrick
      │        │
      │        ├─ Implicit?
      │        │     ├─ near-incompressible / general default ─► -formulation bbar
      │        │     ├─ coarse bending / thin ───────────────► -formulation ssp  (or eas)
      │        │     └─ low-ν, no bending, want 8 clean GPs ──► -formulation std
      │        │
      │        └─ Explicit?
      │              ├─ want 8 damage points ──► -formulation uri -hourglass physical -lumped
      │              └─ cheapest (1 mat eval) ─► -formulation ssp -lumped
      │
      └─ Geometry NOT hex-meshable (organic/CAD) ───► BezierTet10
               ├─ near-incompressible ──► -bbar
               ├─ explicit ─────────────► default lumped mass (no flag needed)
               └─ implicit accuracy ────► -cMass
      │
      ▼  (then, for BezierTet10 or LadrunoBrick, layer the kinematic regime)
   -geom linear   (default; small strain/rotation)
   -geom corot    (large rotation, small strain)
   -geom finite   (large strain; needs FiniteStrainNDMaterial; -bbar ⇒ F-bar ⇒ unsymmetric solver)
```

> [!example] The procedure in one sentence
> *2D → `BezierTri6`. 3D & hex-meshable → `LadrunoBrick` (`bbar` implicit /
> `uri+physical` or `ssp` explicit). 3D & not hex-meshable → `BezierTet10`. Then
> add `-bbar` for incompressibility and `-geom corot|finite` for large
> deformation.*

---

## 4 · Intended-use profiles (per element)

### 4.1 `BezierTri6` — quadratic 2D continuum
> [!info] Reach for it when
> - The model is **2D** (plane stress / plane strain / a representative slice).
> - You want **smooth quadratic stress fields** without the parasitic shear of a
>   bilinear quad, or you are running **explicit 2D dynamics** and need a
>   non-negative lumped mass (its defining feature, Kadapa 2018 Eq. 56).
> - You already have a **quadratic-Lagrange (T6) mesh** from gmsh/apeGmsh — the
>   element consumes it directly via the per-edge control-point map.

**Command**
```python
element('BezierTri6', tag, n1,n2,n3,n4,n5,n6, thick, type, matTag
        [, '-bbar']            # plane strain only (rejected for PlaneStress)
        [, '-cMass']           # consistent mass; default is lumped ρVe/6
        [, '-rho', r]          # else from material
        [, '-pressure', p]     # volume hack (legacy)
        [, '-bodyForce', b1, b2])
```
- Node order: corner1, corner2, corner3, mid12, mid23, mid31.
- `type` ∈ {`PlaneStress`, `PlaneStrain`}. **`-bbar` + `PlaneStress` is a no-op**
  (warned + ignored): plane stress has no volumetric constraint.
- **No** `-geom` — small-displacement only. For large 2D deformation there is no
  Bézier path yet; use a different element or extrude to 3D + `BezierTet10`.
- Straight-sided only (mid-edge nodes at edge midpoints); fixed 3-pt quadrature.

### 4.2 `BezierTet10` — quadratic 3D tetrahedron
> [!info] Reach for it when
> - The geometry is **3D and hard to hex-mesh** (CAD imports, organic shapes,
>   re-entrant features) — automatic tet meshing + quadratic accuracy.
> - **3D explicit dynamics** with quadratic tets (the negative-lumped-mass
>   problem of standard Tet10 is exactly what the Bernstein basis fixes, Eq. 57).
> - **Near-incompressible 3D** (thick-walled sphere is the canonical case): the
>   3D B-bar (Eq. 45) stays optimal as `ν → 0.5` where linear tets lock badly.

**Command**
```python
element('BezierTet10', tag, n1..n10, matTag
        [, '-bbar']                              # near-incompressibility
        [, '-cMass']                             # consistent; default lumped ρVe/10
        [, '-rho', r]
        [, '-bodyForce', b1, b2, b3]
        [, '-pressure', p]                       # volume hack; not under corot/finite
        [, '-geom', <linear|corot|finite>]       # default linear
        [, '-fbar', <centroid|mean_dilatation>]) # only with -bbar -geom finite
```
- Node order (TenNodeTetrahedron): v1 v2 v3 v4, e12 e23 e13 e14 e34 e24.
- `-geom corot` = EICR large-rotation/small-strain (see [[14_bezier_tet10_corot]]);
  `-geom finite` = updated-Lagrangian large strain (needs a `FiniteStrainNDMaterial`).
- `-geom finite -bbar` = **F-bar** (`-fbar centroid` is the dSNPO §15.1 form;
  `mean_dilatation` is the volume-averaged variant) → **unsymmetric tangent**.
- `-pressure` is rejected under `corot`/`finite`. Straight-sided only; 4-pt rule.

### 4.3 `LadrunoBrick` — unified 3D hexahedron
> [!info] Reach for it when
> - The geometry **hex-meshes cleanly** and you want the best accuracy-per-DOF.
> - You want **one element** to span full integration, B-bar, reduced
>   integration + hourglass control, SSP, and enhanced assumed strain — switched
>   by a flag, not by rebuilding the model with a different class.
> - You need a **cheap explicit hex** (`uri`/`ssp`) or **large-deformation hex**
>   (`corot`/`finite` + F-bar).

**Command** (full menu — see [[LadrunoBrick_reference]] for the deep treatment)
```python
element('LadrunoBrick', tag, n1..n8, matTag
        [, '-formulation', <std|bbar|uri|ssp|eas>]      # default std
        [, '-geom',        <linear|corot|finite>]       # default linear
        [, '-hourglass',   <stiffness|physical|viscous> [, coeff]]  # uri only
        [, '-lumped'] [, '-b', bx,by,bz] [, '-damp', dampTag])
```

| `-formulation` | What it cures | Use it for |
|---|---|---|
| `std` | nothing (full 2×2×2) | low-ν, no bending; 8 clean damage points; correctness anchor (== `Brick`) |
| `bbar` | volumetric locking | **default implicit**; near-incompressible (== `bbarBrick`) |
| `uri` | shear + volumetric (1-pt) | **explicit** hex (`+physical` for 8 damage pts, `+stiffness`/`viscous`) |
| `ssp` | shear + volumetric, all ν | coarse bending, `ν→0.5`, cheapest stable choice (== `SSPbrick`) |
| `eas` | parasitic-bending locking | bending-dominated implicit (true Simo–Rifai, ADR 19; small strain) |

- `corot`/`finite` support **`std`/`bbar` only** (uri/ssp/eas under large
  deformation are reserved). `-damp` is honoured by `std`/`bbar` only.
- `uri -hourglass viscous` is **explicit-only** (rate damping, rank-deficient
  tangent under an implicit/eigen solver).

---

## 5 · Cross-cutting modeling guidance

### 5.1 Mass & explicit dynamics
- The Bézier elements default to **lumped** mass (positive by construction) — do
  **not** pass a flag for explicit; pass `-cMass` only when you want consistent
  mass for an implicit modal/accuracy run.
- `LadrunoBrick` defaults to **consistent** mass; pass **`-lumped`** for explicit.
  Its trilinear row-sum lumping is non-negative (unlike serendipity hexes).
- Drive `dt = 0.9 · criticalTimeStep()`. Damping shrinks `dt_cr` in explicit mode
  — see [[central_difference_ladruno_guide]] / [[12_damping_channels]].

### 5.2 Near-incompressibility — pick the right cure
| Element | Cure | Note |
|---|---|---|
| `BezierTri6` | `-bbar` | plane strain only |
| `BezierTet10` | `-bbar` (+ `-fbar` for finite) | |
| `LadrunoBrick` | `-formulation bbar` or `ssp` | `ssp` is well-conditioned across all ν |

Full integration of a near-incompressible material **without** a cure locks —
displacements come out far below the analytic answer. This is the single most
common modeling mistake (Axis E).

### 5.3 Large deformation
- `corot` (large rotation, small strain) and `finite` (large strain) are
  available on `BezierTet10` and `LadrunoBrick`. **`BezierTri6` has neither.**
- `finite` requires a `FiniteStrainNDMaterial`; a small-strain material under
  `-geom finite` is rejected at parse time, and a finite-strain material under
  `-geom linear|corot` is **also** rejected (it is driven only by `setTrialF`).
- F-bar (`-geom finite -bbar`) ⇒ unsymmetric solver. Combined-hardening
  plasticity is non-objective under large rotation (dSNPO §14.11) — a known
  limit, pinned as xfail for `LadrunoJ2`-finite.

### 5.4 Softening / concrete
- All three elements feed a geometry-true `lch` to the material's crack band
  (`sqrt(2A)` / `cbrt(6V)` / `cbrt(V)`), so dissipation is **mesh-objective**.
- Brick + `ASDConcrete3D`: prefer `std`/`bbar`; `ssp`/`uri` work (damage-scaled
  `Kstab`) but **watch `hourglassEnergy`**. Full recipe in
  [[11_brick_asdconcrete_integration]].

### 5.5 Embedded reinforcement
All three elements can act as the **host** for `LadrunoEmbeddedRebar` (the rebar
node ties to the host via host shape functions; pass `-host hostEleTag` and the
host's node list *is* the embedding cage). See
[[20_ladruno_embedded_reinforcement_adr]].

### 5.6 Recording
All three implement the Ladruno recorder contract (self-declared basis/topology)
and the standard `stresses`/`strains`/`material N` response tree, so they record
zero-edit through the Ladruno/MPCO recorders.

---

## 6 · Quick-start recipes

```python
# ── 2D plane-strain, near-incompressible, implicit ──────────────
ops.element('BezierTri6', 11, *n6, 1.0, 'PlaneStrain', mat, '-bbar', '-cMass')

# ── 3D organic geometry, implicit static ────────────────────────
ops.element('BezierTet10', 21, *n10, mat)                 # add -bbar if ν→0.5

# ── 3D explicit impact on a tet mesh ────────────────────────────
ops.element('BezierTet10', 22, *n10, mat)                 # lumped mass is default
ops.integrator('CentralDifferenceLadruno')
dt = 0.9 * ops.criticalTimeStep()

# ── 3D hex, implicit push-over (recommended baseline) ───────────
ops.element('LadrunoBrick', 31, *n8, mat, '-formulation', 'bbar')

# ── 3D hex, coarse bending / ν→0.5 ──────────────────────────────
ops.element('LadrunoBrick', 32, *n8, mat, '-formulation', 'ssp')

# ── 3D hex, explicit dynamic ────────────────────────────────────
ops.element('LadrunoBrick', 41, *n8, mat,
            '-formulation', 'uri', '-hourglass', 'physical', '-lumped')

# ── 3D large strain (hex or tet) with a finite-strain material ──
ops.nDMaterial('LogStrain', 9, baseJ2Tag)
ops.element('LadrunoBrick', 51, *n8, 9, '-formulation', 'bbar', '-geom', 'finite')
ops.system('FullGeneral')        # F-bar tangent is unsymmetric
```

---

## 7 · Anti-patterns & common mistakes

> [!warning] Watch for these
> - **Near-incompressible material on a fully-integrated element** (`std` brick,
>   plain Tet10, `BezierTri6` without `-bbar`) → volumetric locking. Add the cure
>   (Axis E / §5.2).
> - **`-bbar` to fix bending** on the brick → wrong knob. B-bar cures volumetric
>   locking only; use `ssp`/`eas` (or a quadratic Bézier element) for bending.
> - **Explicit run with `LadrunoBrick` but no `-lumped`** → consistent mass, no
>   cheap diagonal solve. Conversely, passing `-lumped`/forgetting `-cMass` is
>   irrelevant for the Bézier elements (lumped is already their default).
> - **`-geom finite` with a small-strain material** (or a `FiniteStrainNDMaterial`
>   under `linear`/`corot`) → rejected at parse time. Match the material to the
>   geometry method.
> - **F-bar with a symmetric solver** (`BandSPD`/`ProfileSPD`) → the unsymmetric
>   coupling term is dropped and Newton may not converge. Use `FullGeneral`/
>   `UmfPack`/`SparseGEN`.
> - **Reaching for `BezierTet10` as the default 3D element** → if the geometry
>   hex-meshes, `LadrunoBrick` is more accurate per DOF. Tets are for geometry
>   that forces them (Axis B).
> - **Large-deformation 2D with `BezierTri6`** → it has no `corot`/`finite`. Move
>   to 3D + `BezierTet10`, or another element.

---

## 8 · References & maintenance

**Per-element references:** [[LadrunoBrick_reference]] ·
[[04_bezier_elements|BezierTri6 ADR]] · [[06_bezier_tet10|BezierTet10 ADR]] ·
[[14_bezier_tet10_corot|BezierTet10 corot]].

**Companion notes:** [[11_brick_asdconcrete_integration]] (concrete/softening) ·
[[solid_transformation_wrapper]] (linear/corot/finite seams) ·
[[09_finite_strain_material_wrapper]] (LogStrain / `setTrialF`) ·
[[20_ladruno_embedded_reinforcement_adr]] (embedded rebar host) ·
[[central_difference_ladruno_guide]] / [[explicit_bathe_integrators_guide]] (explicit) ·
[[12_damping_channels]] · [[LEDGER_implementations]] (authoritative tags/PRs) ·
[[LEDGER_quirks]].

**Source:** `SRC/element/bezierTriangle/` (`33000`) ·
`SRC/element/bezierTetrahedron/` (`33001`) · `SRC/element/ladrunoBrick/` (`33002`)
· `SRC/classTags.h:918-920`.

**Theory:** Kadapa, C. (2018) *Novel quadratic Bézier triangular and tetrahedral
elements …*, IJNME 117:543 (Bézier elements). Hughes (1980) — B-bar. de Souza
Neto, Perić & Owen — *Computational Methods for Plasticity* §15.1 (F-bar), §14
(finite-strain). Simo & Rifai (1990) — EAS. Flanagan & Belytschko (1981) —
hourglass control. (Full per-element bibliographies in the deep docs.)

---

> [!info] Maintenance
> This is a **living selection guide**. When an element gains/loses a formulation,
> geometry method, or constraint, update the §1 table and the affected profile in
> §4 here, **and** the per-element deep doc, [[LEDGER_implementations]], and the
> banner (`Ladruno_scripts/banner_features.txt`). This guide stays at the
> *selection* altitude — push new theory/implementation detail down into the
> per-element references.
