---
title: Lysmer–Kuhlemeyer & ASD free-field absorbing boundaries — deep guide
project: Ladruno
status: reference
owner: nmora
tags:
  - ssi
  - absorbing-boundaries
  - wave-propagation
  - apegmsh
  - reference
---

# Lysmer–Kuhlemeyer & ASD free-field absorbing boundaries — deep guide

A **source-accurate, formulation-level** study of OpenSees' viscous-dashpot family
of wave-truncation boundaries — `LysmerTriangle`, `ASDAbsorbingBoundary2D`, and
`ASDAbsorbingBoundary3D` — written to three audiences at once:

1. a researcher who needs to understand **what the element actually computes**
   (theory + the assembled M/C/K and resisting force, term by term);
2. a developer who needs the **OpenSees architecture and implementation** (the
   stage machine, the node sorting, the static condensation, the penalty
   bookkeeping);
3. the **apeGmsh team**, who must **emit and drive** these elements from a mesh —
   §9 is a self-contained generator contract (ghost-layer construction, `btype`
   strings, staging, input convention, what to read back).

All `file:line` citations are against this worktree's `SRC/`. This is a **read /
document** pass — no source was modified. The broader survey of *every*
wave-truncation facility (PML, DRM) lives in
[[absorbing_boundaries_and_pml_guide]]; **this** guide is the deep dive on the
dashpot / free-field family and is the normative reference for an apeGmsh
implementation.

> **Upstream, not fork-authored.** These three elements are stock OpenSees
> (`LysmerTriangle` by J. A. Abell, UANDES; the two `ASDAbsorbingBoundary*` by
> M. Petracca & G. Camata, ASDEA). Nothing here is Ladruno-only — but the SSI /
> explicit-dynamics work on this fork leans on them, and apeGmsh must generate
> them, so they get the same deep treatment as a fork feature.

---

## 1. TL;DR — what these elements really are

The single most important thing to understand: **the two ASD elements are not bare
Lysmer dashpots.** They implement the **Lysmer free-field boundary with a
compliant base** — the same architecture FLAC and PLAXIS use for SSI — where each
element straddles the truncation interface and carries *two* columns of nodes:

- **SS (soil-side) nodes** — shared with your actual soil mesh (the real boundary
  nodes of the truncated domain).
- **FF (free-field) nodes** — extra "ghost" nodes *outside* the mesh that run an
  independent **1-D free-field column** (a built-in site-response solver).

In the absorbing stage the element superposes **four** mechanisms:

| # | Mechanism | Lives on | Code (2D) |
|---|---|---|---|
| 1 | **Free-field column** (1-D site response: K, M, Rayleigh C) | FF nodes | `addKff`/`addMff`/`addCff` |
| 2 | **Free-field → soil traction** `σ·n` (injects the free-field motion) | SS nodes | `addKffToSoil`/`addRffToSoil` |
| 3 | **Lysmer–Kuhlemeyer dashpots on the *relative* velocity** `V_ff−V_ss` | FF↔SS | `addClk`/`addRlk` |
| 4 | **Compliant base** input `2ρV·v_incident` (bottom faces only) | SS nodes | `addBaseActions` |

Mechanism **3 is the key idea**: because the dashpot acts on the *relative*
velocity, the free-field component cancels and **only the scattered (outgoing)
field is absorbed** — the free-field passes through the boundary unattenuated. A
bare Lysmer dashpot (mechanism 3 against a *fixed* reference) cannot do that; it
eats the free-field too, which is why a bare dashpot needs the input re-injected as
a force and is exact only at normal incidence.

`LysmerTriangle` *is* the bare dashpot (mechanism 3 against ground, + an optional
gravity-hold spring). It is lighter, triangular, and the clean choice for explicit
dynamics; the ASD elements are the production free-field boundary.

**Mental model for a complete model:**

```
Interior (structure + nonlinear near-field soil)
  └─ [optional DRM ring — injects a full 3-D regional field]
       └─ ASDAbsorbingBoundary ghost ring:
            • sides  (L/R/F/K): free-field column + relative-velocity dashpots
            • bottom (B):       compliant base + 2ρV·v_incident input (VELOCITY)
            • edges/corners:    ONE element, combined btype ("BLF"…)
       └─ Rayleigh damping assigned to the absorbing elements (NOT optional)
Workflow: build(stage 0) → gravity static → loadConst → setParameter stage=1 → transient
```

---

## 2. Theory

### 2.1 The Lysmer–Kuhlemeyer viscous boundary (1969)

A plane body wave normally incident on a truncation plane is perfectly absorbed if
the boundary applies a traction equal to the medium impedance times the particle
velocity:

```
t_n = ρ · V_p · v_n        (normal / P-wave),   V_p = √((λ+2μ)/ρ)
t_s = ρ · V_s · v_t        (tangent / S-wave),  V_s = √(μ/ρ)
```

Integrated over a boundary facet of tributary area `A`, this is a set of nodal
dashpots `c_n = ρV_p·A`, `c_t = ρV_s·A`. The absorption is **exact only at normal
incidence**; at incidence angle θ a fraction of energy reflects, and the
performance degrades for surface waves and at low frequency (wavelength ≫ domain).
This is the first-order absorbing boundary — cheap, explicit-friendly, but leaky.

### 2.2 Why a free-field column (the lateral-boundary problem)

For an SSI box, the lateral boundaries must do two contradictory things: **absorb**
the scattered field radiating away from the structure, **and reproduce** the
far-field 1-D site response (the soil would be shaking even with no structure). A
bare dashpot can only do the first.

The **free-field boundary** solves this by attaching, at each lateral boundary
element, a 1-D soil column (the "free-field column") that is driven by the same base
motion and develops the correct site response. Then:

- the column's stress is applied as a traction on the soil-side nodes (drive the
  domain with the free-field), and
- dashpots connect the soil nodes to the column nodes and act on the **velocity
  difference**, absorbing only the part of the motion that deviates from the
  free-field — i.e. the scattered field.

This is the Lysmer free-field boundary (Zienkiewicz et al.; widely known from the
FLAC/PLAXIS manuals; see also Nielsen 2006, Kontoe et al. 2009). The ASD elements
are a finite-element packaging of exactly this idea, with the free-field column
*embedded in the boundary element itself* rather than as a separate attached mesh.

### 2.3 The compliant base (bottom boundary)

At the base you usually want to **input** a motion *and* absorb the down-going
reflected field. The **compliant base** (a.k.a. Joyner; "within" motion) does both:
apply, at the base nodes, a force equal to `2 · ρV · v_incident` while a Lysmer
dashpot sits in parallel. The dashpot absorbs the reflected half; the factor of 2
accounts for it, so the *resulting* base motion equals the incident ("within")
motion. The corollary: **you feed the base the incident (outcrop ÷ 2) velocity
history**, not the full outcrop and not an acceleration.

### 2.4 Static stage (gravity)

None of the above is correct under gravity — you do not want the boundary radiating
or free-fielding while you establish the initial stress state. So the ASD elements
run a **two-stage** life: a **static-constraint stage** (the boundary behaves as a
roller / tied boundary, holding the domain for gravity) and then a one-way switch
to the **absorbing stage** for the dynamic run, carefully transferring the gravity
reactions so the static stress state survives the switch (§5.2, §6.2).

---

## 3. Class tags, files, command registration

| Element | Class tag | Files | Author |
|---|---|---|---|
| `LysmerTriangle` | `ELE_TAG_LysmerTriangle = 185` | [LysmerTriangle.{h,cpp}](SRC/element/absorbentBoundaries/LysmerTriangle.cpp) | J. A. Abell (UANDES) |
| `ASDAbsorbingBoundary2D` | `ELE_TAG_ASDAbsorbingBoundary2D = 219` | [ASDAbsorbingBoundary2D.{h,cpp}](SRC/element/absorbentBoundaries/ASDAbsorbingBoundary2D.cpp) | M. Petracca, G. Camata (ASDEA) |
| `ASDAbsorbingBoundary3D` | `ELE_TAG_ASDAbsorbingBoundary3D = 220` | [ASDAbsorbingBoundary3D.{h,cpp}](SRC/element/absorbentBoundaries/ASDAbsorbingBoundary3D.cpp) | M. Petracca, G. Camata (ASDEA) |
| `LysmerVelocityLoader` | (elemental load) | [LysmerVelocityLoader.{h,cpp}](SRC/domain/load/LysmerVelocityLoader.cpp) | — |

- Build: [SRC/element/absorbentBoundaries/CMakeLists.txt](SRC/element/absorbentBoundaries/CMakeLists.txt).
- Modern command registration (openseespy/Tcl `element`): the `OPS_*` factory
  functions at the top of each `.cpp` (`OPS_LysmerTriangle`,
  `OPS_ASDAbsorbingBoundary2D`, `OPS_ASDAbsorbingBoundary3D`).
- Broker (parallel/clone): `FEM_ObjectBrokerAllClasses.cpp`.

These are all **upstream** tags (<33000), not fork-allocated.

---

## 4. `LysmerTriangle` — the bare dashpot facet

3-node triangular facet, **3 DOF/node = 9 DOF**, placed as a skin on a truncation
surface. Provides only a **damping** matrix (and an optional gravity-hold spring);
**no mass** in any stage (inertia comes from the adjacent soil solids).

### 4.1 Command / API

```tcl
element LysmerTriangle $eleTag $iNode $jNode $kNode  $rho $Vp $Vs  <$length> <$stage>
```
([parser, LysmerTriangle.cpp:67](SRC/element/absorbentBoundaries/LysmerTriangle.cpp))

| arg | meaning |
|---|---|
| `$rho` | mass density of the adjacent soil |
| `$Vp`, `$Vs` | P- and S-wave velocities of the adjacent soil |
| `$length` *(opt)* | explicit length scale for the gravity-hold spring; if `0`/omitted the element uses the average side length |
| `$stage` *(opt, default 0)* | operating mode, **set at construction** (not switchable live) |

### 4.2 Geometry (the local triad)

[`UpdateBase`, line 264](SRC/element/absorbentBoundaries/LysmerTriangle.cpp) builds
an orthonormal triad from the facet edges `g1 = x_j−x_i`, `g2 = x_k−x_i`:

```
n̂ = g1 × g2          (facet normal),    A = ‖g1 × g2‖ / 2   (facet area)
t̂ = g1               (tangent-1)
ŝ = n̂ × g1           (tangent-2)
```
The constant `Bmat` (9×3 of 0.5's, set in `setDomain`) lumps the facet quantities
equally to the three nodes.

### 4.3 Damping matrix (the dashpot) — stages 0/2/3

[`getDamp`, line 376](SRC/element/absorbentBoundaries/LysmerTriangle.cpp):

```
C_local = diag(ρV_s, ρV_s, ρV_p)            // in (t, s, n)
C_global = Bᵀ · ( Tᵀ C_local T · A ) · B    // rotate to global, lump to nodes
```
with `T` the (t̂, ŝ, n̂) rotation. Each node receives the Lysmer dashpot scaled by
the tributary area. This is the standard viscous boundary.

### 4.4 Gravity-hold spring — stages 1/2

[`getTangentStiff`, line 304](SRC/element/absorbentBoundaries/LysmerTriangle.cpp):
a **normal-direction-only** elastic spring (the two shear terms `K(0,0)`/`K(1,1)`
are deliberately commented out):

```
G = ρV_s²,   M = ρV_p²,   E = G(3M − 4G)/(M − G),   K_local(n,n) = E/L
```
`E` is just Young's modulus recovered from `G` and the constrained modulus `M`
(since `λ = M − 2G`, `E = μ(3λ+2μ)/(λ+μ)`). The spring exists only to restrain the
boundary in the normal direction during a static/gravity step.

### 4.5 Stage semantics

| `stage` | behaviour |
|---|---|
| **0** (default) | **pure dashpot** (damping only) — the dynamic absorbing mode |
| 1 | **pure normal spring** (gravity hold), no damping |
| 2 | spring **and** dashpot |
| 3 | dashpot + pass the committed spring force straight through as a reaction |

Unlike the ASD elements, `LysmerTriangle` has **no live stage switch** — you choose
the stage at construction. The intended workflow is to build the gravity model with
stage-1 facets, then rebuild/recreate them as stage-0 for the dynamic run (or use a
separate gravity restraint).

### 4.6 Velocity input & a residual subtlety

The companion **`LysmerVelocityLoader`** (an `ElementalLoad`, ctor `(tag, eleTag,
dir)`) feeds a prescribed boundary velocity into the element in direction `dir`
(0/1/2), accumulated into `gnd_velocity` via `addLoad` and reset each `zeroLoad`.
This is the "dashpot + incident velocity" input path without full DRM.

**Subtlety worth pinning** ([`getResistingForceIncInertia`, line 500](SRC/element/absorbentBoundaries/LysmerTriangle.cpp)):
the residual is formed with `velocity = 0*v_node + gnd_velocity` — it **literally
zeroes the node velocity** and keeps only the prescribed ground velocity. So
`getResistingForceIncInertia` returns only the *incident-wave forcing* `C·v_gnd`;
the reactive `C·v_node` dashpot force is delivered through the **assembled global
damping matrix** (`getDamp()`), which the transient integrator multiplies by the
nodal velocity. The exact residual-vs-tangent split is integrator-dependent — if it
ever matters for a specific explicit scheme, confirm with a single-element
velocity-step test rather than trusting the comment.

---

## 5. `ASDAbsorbingBoundary2D` — the production free-field boundary (2D)

A **4-node quad** straddling the interface, **2 DOF/node = 8 DOF**. Plane-strain
soil with shear modulus `G`, Poisson `v`, density `rho`, out-of-plane `thickness`.

### 5.1 Command / API

```tcl
element ASDAbsorbingBoundary2D $tag $n1 $n2 $n3 $n4  $G $v $rho $thickness $btype  <-fx $tsxTag> <-fy $tsyTag>
```
([parser, ASDAbsorbingBoundary2D.cpp:142](SRC/element/absorbentBoundaries/ASDAbsorbingBoundary2D.cpp))

| arg | meaning |
|---|---|
| `$n1..$n4` | the 4 nodes — **order does not matter**, the element sorts them by coordinate (§5.3) |
| `$G $v $rho` | shear modulus, Poisson's ratio, density of the adjacent soil |
| `$thickness` | out-of-plane thickness |
| `$btype` | boundary-position string: `"B"` bottom, `"L"` left, `"R"` right — **combinable** for corners (`"BL"`, `"BR"`) |
| `-fx/-fy $ts` *(opt, bottom only)* | `TimeSeries` providing the base **velocity** action per direction (compliant-base input) |

The `-fx/-fy` options are parsed **only when `$btype` contains "B"** — base actions
exist only on the bottom.

### 5.2 The stage machine

`enum StageType { Stage_StaticConstraint = 0, Stage_Absorbing = 1 }`
([ASDAbsorbingBoundary2D.h:40](SRC/element/absorbentBoundaries/ASDAbsorbingBoundary2D.h)).

Live, one-way switch via a parameter named `"stage"`:
```tcl
setParameter -val 1 -ele $tag1 $tag2 ... stage    ;# 0 -> 1 only
```
[`updateParameter`, line 816](SRC/element/absorbentBoundaries/ASDAbsorbingBoundary2D.cpp)
rejects anything but `0→1` with `exit(-1)`. The switch runs
[`updateStage`, line 1005](SRC/element/absorbentBoundaries/ASDAbsorbingBoundary2D.cpp):

```
m_R0 = -(stage-0 penalty reactions)   // saved, re-applied in stage 1
m_U0 += current displacements         // datum captured at the switch
m_stage = Stage_Absorbing
```

Two consequences make the gravity→dynamics handoff clean:

- **`m_R0`** is re-applied as an external load in stage 1
  ([`addRReactions`, line 1215](SRC/element/absorbentBoundaries/ASDAbsorbingBoundary2D.cpp)),
  so releasing the constraints does **not** dump the gravity stress.
- **`m_U0`** is the datum: `getDisplacement()` returns `U − U0`. Every free-field /
  dashpot term therefore runs on the **dynamic increment** — gravity displacements
  never pollute the free-field stress.

Query the current stage with `eleResponse $tag stage` (also `G`, `v`, `rho`, `E`;
[`setResponse`, line 861](SRC/element/absorbentBoundaries/ASDAbsorbingBoundary2D.cpp)).

### 5.3 Internal node sort & DOF map

On first `setDomain`, the 4 nodes are sorted ([`SorterLeft`/`SorterRight`, line 98](SRC/element/absorbentBoundaries/ASDAbsorbingBoundary2D.cpp))
left-to-right then bottom-to-top (right-to-left for a `RIGHT` boundary), so the
**local layout is deterministic regardless of input order**:

```
local 0,1 = OUTER (free-field, FF) column   (further from the soil)
local 2,3 = INNER (soil-side,   SS) column   (shared with the soil mesh)
```
`m_dof_map(8)` maps local DOFs → element matrix positions; `getElementSizes`
([line 986](SRC/element/absorbentBoundaries/ASDAbsorbingBoundary2D.cpp)) recovers
`lx` (across-thickness), `ly` (along-boundary), and the outward-normal sign
`nx = +1` (LEFT) / `−1` (RIGHT). **Geometry must be axis-aligned** — sizes are raw
coordinate differences.

### 5.4 Stage 0 — static constraint (gravity)

[`addKPenaltyStage0`, line 1033](SRC/element/absorbentBoundaries/ASDAbsorbingBoundary2D.cpp)
enforces, by **penalty**:

- **SP**: fix the normal DOF at every node (Uy on a bottom face / Ux on a vertical
  face) → roller.
- **MP** (`equalDOF`): tie the tangential DOF of the FF nodes to the SS nodes → the
  ghost column tracks the soil.

Penalties auto-scale to the material
([`penaltyFactor`, line 1019](SRC/element/absorbentBoundaries/ASDAbsorbingBoundary2D.cpp)):
`sp = 10^(oom+8)`, `mp = 10^(oom+3)`, `oom = round(log10(G·t))`. The `+8` keeps `sp`
≈ 1e8× the natural stiffness — strong as a constraint, still inside a double's
~16 significant digits.

### 5.5 Stage 1 — the four mechanisms (vertical face)

**(1) Free-field column** — a 1-D column between the two FF nodes (0,1):
```
λ = 2μv/(1−2v),   μ = G
k_shear  = lx·μ·t/ly          (tangential, S-wave column)   addKff, line 1280
k_axial  = lx·(λ+2μ)·t/ly     (normal,    P-wave column)
M_ff     = ρ·t·lx·ly/2  lumped to the 2 FF nodes            addMff, line 1229
C_ff     = α·M_ff + β·K_ff    (Rayleigh, element's own α,β)  addCff, line 1413
```

**(2) Free-field → soil traction** — the column stress applied on the SS nodes
(non-symmetric, one-way: free-field drives soil)
([`addKffToSoil`, line 1343](SRC/element/absorbentBoundaries/ASDAbsorbingBoundary2D.cpp)):
coefficients `±λ·nx·t/2` (axial) and `±μ·nx·t/2` (shear). This is what injects the
lateral free-field motion into the domain.

**(3) Lysmer dashpots on the relative velocity**
([`getLKcoeff`/`addClk`, line 1472](SRC/element/absorbentBoundaries/ASDAbsorbingBoundary2D.cpp)):
```
ap = ρ·V_p·h·t/2,   as = ρ·V_s·h·t/2,   h = ly      // force = c·(V_ff − V_ss)
```
Acting on `V_ff − V_ss` is the crux — the free-field cancels, only the scattered
field is damped.

### 5.6 Stage 1 — bottom (horizontal) face: a compliant base

On a bottom face the free-field column is **off** (`addKff/addMff/addCff`
early-return on `BND_BOTTOM`). Instead:

- dashpots, with `V_p`/`V_s` **swapped** and `h = lx` (normal is now vertical,
  [line 1490](SRC/element/absorbentBoundaries/ASDAbsorbingBoundary2D.cpp)), and
- **base actions** ([`addBaseActions`, line 1580](SRC/element/absorbentBoundaries/ASDAbsorbingBoundary2D.cpp)):
  ```
  f_x = 2 · v_x · (ρV_s)      // horizontal → S-wave impedance
  f_y = 2 · v_y · (ρV_p)      // vertical   → P-wave impedance
  ```
  `v_x`, `v_y` come from `-fx`/`-fy` `TimeSeries` as **velocity** histories. The
  **factor of 2** is the compliant-base signature (feed the within/incident motion).

### 5.7 Resisting-force assembly (the order of operations)

[`getResistingForceIncInertia`, line 524](SRC/element/absorbentBoundaries/ASDAbsorbingBoundary2D.cpp),
stage 1:
```
R  = addRPenaltyStage1   (bottom-only SP on outer nodes)
   + addRff              (free-field internal force)
   + addRffToSoil        (traction onto soil)
   + addRReactions       (restore stage-0 gravity reactions)
   + addBaseActions      (compliant-base input)
   + addRCff + addRlk    (C·v damping forces: free-field + dashpots)
   + addRMff             (M·a inertia of the free-field column)
```
The tangent (`getTangentStiff`/`getDamp`/`getMass`) returns the matching `K`, `C`,
`M`. **`C·v` appears in both** the residual (here) and the damping matrix (for the
Jacobian) — that is correct, not double-counted: one is the residual, one is the
tangent.

---

## 6. `ASDAbsorbingBoundary3D` — the same idea, generalized

An **8-node hex**, **3 DOF/node = 24 DOF**. Same four mechanisms and same stage
machine; three things make the 3D version more involved.

### 6.1 Command / API

```tcl
element ASDAbsorbingBoundary3D $tag $n1 ... $n8  $G $v $rho $btype  <-fx $tsx> <-fy $tsy> <-fz $tsz>
```
([parser, ASDAbsorbingBoundary3D.cpp:473](SRC/element/absorbentBoundaries/ASDAbsorbingBoundary3D.cpp))

- 8 nodes (order-independent — sorted internally), `G v rho`, **no thickness**
  (3-D), and `$btype` over `"B" "L" "R" "F"(front) "K"(back)`, combinable.
- `-fx/-fy/-fz` base **velocity** actions, bottom faces only.

### 6.2 The 17-way boundary classification

A 3-D absorbing brick can sit on a **side** (1 face), a **vertical edge** (2
vertical faces), a **bottom edge** (bottom + 1 vertical), or a **bottom corner**
(bottom + 2 verticals). `$btype` is bit-flagged into `BND_BOTTOM/LEFT/RIGHT/FRONT/
BACK`; one of five node-sorters is picked per case
([line 742](SRC/element/absorbentBoundaries/ASDAbsorbingBoundary3D.cpp)), and the
Lysmer node-pairs come from `LKselectPairs` with **tributary weights** — `w=2` on
vertical edges, `w=4` on bottom corners, and the free-field mass split `÷2`/`÷4`
([`LK_NODES_*`, line 253](SRC/element/absorbentBoundaries/ASDAbsorbingBoundary3D.cpp)).
**You place ONE element per ghost brick with the combined `btype`** (e.g. `"BLF"`);
the element does the corner/edge tributary bookkeeping internally — never overlap
two absorbing elements on the same brick.

### 6.3 The free-field column is a statically-condensed H8

Where 2D hand-codes the 1-D column, 3D builds the **real 8-node-hex elastic
stiffness** with full 2×2×2 Gauss
([`addKff`, line 1948](SRC/element/absorbentBoundaries/ASDAbsorbingBoundary3D.cpp))
but routes the assembly through
[`ffMapping`, line 1894](SRC/element/absorbentBoundaries/ASDAbsorbingBoundary3D.cpp),
which **slaves the SS face DOFs to the FF face DOFs** (`n2←n0, n3←n1, n6←n4,
n7←n5`). The condensed B-matrix then represents a **pure 1-D free-field column** (no
in-plane gradient). The traction back onto the soil uses a proper Voigt `σ·n`
operator lumped over the face area
([`computeNmatrix`, line 2092](SRC/element/absorbentBoundaries/ASDAbsorbingBoundary3D.cpp)),
which **asserts the boundary normal is ±X or ±Y** (`exit(-1)` otherwise) — vertical
faces must be axis-aligned; the bottom (normal Z) is handled by base actions, not
the N-matrix.

### 6.4 Distortion handling

[`handleDistortion`, line 376](SRC/element/absorbentBoundaries/ASDAbsorbingBoundary3D.cpp)
orthogonalizes the Jacobian columns and rebuilds an **equivalent-volume rectangular
brick** before sorting/sizing, so a moderately skewed mesh still yields clean
`lx,ly,lz` and axis-aligned normals (a singular Jacobian `exit(-1)`s). One latent
oddity: the centroid is computed as `sum/3.0` over 8 nodes
([line 437](SRC/element/absorbentBoundaries/ASDAbsorbingBoundary3D.cpp)) — should be
`/8.0` — but it is **harmless**, because `C` only enters as a common additive offset
that cancels in the size *differences* (`N2−N0`) and in the relative sort order.
Worth knowing it is there.

### 6.5 Base actions (3D)

[`addBaseActions`, line 2543](SRC/element/absorbentBoundaries/ASDAbsorbingBoundary3D.cpp):
```
ap = -ρ·V_p·lx·ly/4,   as = -ρ·V_s·lx·ly/4
f = { 2·as·v_x, 2·as·v_y, 2·ap·v_z }    // X,Y horizontal→S ; Z vertical→P
```
distributed to the top (soil-side) nodes with the side/edge/corner multiplicity
matching §6.2.

---

## 7. Operational architecture — how to actually drive them

### 7.1 The staged analysis (canonical, implicit)

```tcl
model BasicBuilder -ndm 3 -ndf 3
# ... soil solids + ASDAbsorbingBoundary3D ghost ring (created in stage 0) ...
# assign Rayleigh damping to the absorbing elements (see 7.3) — e.g. via a region
region 1 -ele $absorb_tags -rayleigh $alphaM 0.0 $betaK 0.0

# --- gravity (stage 0: boundary held by penalty constraints) ---
constraints Transformation
numberer RCM
system UmfPack
integrator LoadControl 1.0
algorithm Newton
analysis Static
analyze 1
loadConst -time 0.0

# --- flip every absorbing element to absorbing mode (one-way 0->1) ---
setParameter -val 1 -ele $absorb_tags stage

# --- dynamics (stage 1: free-field + dashpots + compliant base) ---
integrator Newmark 0.5 0.25
analysis Transient
analyze $nsteps $dt
```

### 7.2 The input convention (read this twice)

- **Input is a *velocity* time history**, supplied as a `TimeSeries` on `-fx/-fy/
  -fz` of the **bottom** absorbing elements (it is multiplied by `ρV` to form a
  stress).
- It is the **within / incident (outcrop ÷ 2)** motion — the factor of 2 in
  `addBaseActions` is what makes the resulting base motion equal the incident
  motion.
- Sides take **no** input; they free-field automatically from the base motion.
- Alternatively, drive the box with a **DRM ring** just inside the absorbing layer
  (full 3-D incident field) and use the absorbing layer purely to absorb — then you
  do **not** use `-fx/-fy/-fz`. The two input paths are alternatives, not stacked.

### 7.3 Damping — the free-field column needs its OWN Rayleigh

`addCff` reads the **element's own** `alphaM/betaK` (the base-class Rayleigh factors
set by `region ... -rayleigh` or `rayleigh` on the element;
[`getDampParam`, line 1402](SRC/element/absorbentBoundaries/ASDAbsorbingBoundary2D.cpp)
picks `betaK`, else `betaK0`, else `betaKc`). **If you do not assign Rayleigh
damping to the absorbing elements, the free-field column is undamped** and will
resonate as a 1-D elastic column at a different damping than your interior soil →
spurious boundary response. Match the absorbing elements' Rayleigh to the interior
soil's.

### 7.4 Solver notes

- The ASD tangent is **non-symmetric** (the one-way FF→soil coupling). A symmetric
  SPD solver (`BandSPD`/`ProfileSPD`) still runs because the asymmetry is small and
  off the absorbed diagonal, but the honest choice on a stiff/awkward case is
  `system UmfPack`/`FullGeneral`.
- It is a **penalty** formulation — expect a tiny constraint residual; do not stack
  conflicting `fix`/`equalDOF` on the boundary DOFs.
- **PML, by contrast, is implicit-only and needs extra nodal DOFs** — the ASD/Lysmer
  family does **not** (standard ndf 2/3). For explicit dynamics the dashpot family
  is the only option (see [[absorbing_boundaries_and_pml_guide]] §6).

### 7.5 Explicit dynamics

`LysmerTriangle` stage-0 supplies a **diagonal** dashpot (`diag(ρV_s,ρV_s,ρV_p)`)
that folds cleanly into central difference — the structurally clean explicit choice.
`ASDAbsorbingBoundary*` supplies coupling (non-diagonal) `M`/`C`/`K` and a
non-symmetric tangent at the boundary ring — workable explicit but not diagonal. For
the fork's explicit integrators see [[project_robust_central_difference]]; remember
**damping shrinks `dt_cr`** in explicit schemes.

---

## 8. Gotchas worth pinning

1. **The free-field column needs its own Rayleigh damping** (§7.3). The single most
   overlooked detail — without it the boundary free-field is undamped.
2. **Input is *velocity*, and the within/÷2 motion** (§7.2). Feeding acceleration or
   full outcrop is wrong by a derivative and a factor of 2.
3. **One element per ghost brick, combined `btype` at edges/corners** (§6.2). The
   element does tributary weighting; overlapping elements double-count.
4. **Boundaries must be axis-aligned.** Sizes are coordinate differences; 3D
   `computeNmatrix` `exit(-1)`s if a vertical normal isn't ±X/±Y.
5. **The stage switch is live and one-way** (`0→1`). Build in stage 0, gravity, then
   `setParameter -val 1 ... stage`. Any other transition `exit(-1)`s.
6. **`btype` orients the element** — get `"L"/"R"/"B"/"F"/"K"` (and the corner
   combination) wrong and the dashpots/free-field silently mis-orient.
7. **Non-symmetric, penalty tangent** (§7.4) — solver-aware.
8. **`LysmerTriangle` has no live stage switch** and **no mass** — gravity vs.
   dynamic is a construction-time choice, and inertia must come from the soil.
9. **Latent `/3.0` centroid in 3D `handleDistortion`** (§6.4) — harmless (cancels in
   differences), but don't "fix" it without re-deriving; it isn't used as a true
   centroid.
10. **The boundary is LINEAR ELASTIC — keep nonlinearity in the interior** (§9.6).
    The free-field column and dashpots use the constant `G,v,rho` you pass and never
    yield. Place the truncation far enough out that the boundary soil stays
    low-strain; a strongly-straining far field makes the (linear) free-field diverge
    from the interior → spurious reflections. Pass the **small-strain `G₀`**, not a
    secant `G`; for layered/pressure-dependent soil use **per-element / depth-graded**
    `G₀` (§9.6). This is the genuine limitation vs FLAC/PLAXIS nonlinear free-fields.
11. **It is an *element*, not a material on a soil brick, and there is no tie** (§9.0).
    `G,v,rho` are numbers in the element command (no `nDMaterial`); coupling is by
    **shared nodes** only (no `equalDOF`/host reference). No soil solids live in the
    ghost layer.

---

## 9. apeGmsh implementation contract — generating & driving the boundary

This section is the deliverable for the apeGmsh team. The irreducible insight:

> **An ASD absorbing boundary is a one-element-thick *ghost layer* wrapped around
> the soil domain.** apeGmsh already owns the geometry and the boundary mesh, so it
> is the natural place to extrude that ghost layer, set each element's `btype` from
> which face(s) it touches, and emit the elements + the staging workflow.

### 9.0 The mental model — what it IS and what it is NOT

These four points are exactly where a first-time implementer (or user) goes wrong;
the generator's design and its docstrings should make them impossible to get wrong.

- **It IS a standalone finite element, NOT a material on a soil brick.**
  `ASDAbsorbingBoundary2D/3D` is a quad/hex element with its own M/C/K. There is **no
  absorbing `nDMaterial`** — the soil properties `G, v, rho` are passed as **plain
  numbers in the element command**. The generator must **not** create a soil solid +
  attach a material in the ghost layer; it creates the ASD element directly. There
  are **no soil solid elements** in the ghost band.
- **There is NO tie / NO host reference.** Unlike `ASDEmbeddedNodeElement` or
  `LadrunoEmbeddedRebar` (which store a host-element tag + weights), the absorbing
  element couples to the soil **only through shared nodes** — its inner face *is* the
  soil boundary nodes. The generator emits **no `equalDOF`, no tie element, no
  constraint**. (The element runs its *own* internal penalty constraints in stage 0,
  but those are inside the C++; nobody authors them.)
- **One element per ghost *cell*, not per face.** A face is a *grid* of ghost cells;
  each cell → one ASD element. A 10×10 base ⇒ 100 `"B"` elements, not one "bottom
  slab." Think "tile the shell with elements," not "5 big slabs."
- **You must fill the edges and corners** (the "diagonal" seams between face patches)
  with their own ghost cells, each carrying the **combined** `btype` — see §9.3. The
  taxonomy for a box open at the top is **17** distinct `btype` strings (5 faces +
  8 edges + 4 bottom corners); the top is a free surface with **no** absorbing layer.

### 9.1 What apeGmsh must generate

For a soil domain meshed with solids (quads in 2D, hexes/`SSPbrick`/`stdBrick` in
3D), to wrap face set(s) `F`:

1. **Extrude a ghost node layer** one element-width *outward* from each boundary
   face. The **inner** face of the ghost layer **reuses the existing soil boundary
   nodes** (shared); the **outer** face is **new ghost nodes** that belong only to
   the absorbing elements (not connected to any soil element). Make the ghost layer
   **conforming** (adjacent ghost bricks share their common nodes, exactly like the
   soil mesh) and **match the adjacent soil element size** (the ghost width `lx`
   sets the free-field column stiffness `∝ lx/ly`).
2. **One absorbing element per ghost brick** (`ASDAbsorbingBoundary2D` per ghost
   quad / `ASDAbsorbingBoundary3D` per ghost hex). **Node order is free** — the
   element sorts internally — but all 4/8 corners of the ghost brick must be passed.
3. **`btype` from face membership** (see 9.3). Edges/corners get the **combined**
   string and a **single** element (do not place two).
4. **Properties `G, v, rho`** passed as **numbers** (no `nDMaterial` object), equal
   to the adjacent soil's **small-strain** values (convert from `E, v` if needed:
   `G = E/(2(1+v))`). 2D also needs `thickness` (match the soil's plane-strain
   thickness). For layered / depth-graded / nonlinear soil these are **per-element**,
   not one global value — see **§9.6**.
5. **Node ndf** ≥ 2 (2D) / ≥ 3 (3D) on the ghost nodes — plain displacement DOFs.
6. **Rayleigh damping** assigned to the absorbing elements (emit a `region ... -ele
   <absorb> -rayleigh ...` matching the interior soil) — **required**, not optional.

### 9.2 Geometry rules apeGmsh must honor

- **Axis-aligned faces only.** Vertical absorbing faces must be ⟂ X or ⟂ Y; the base
  ⟂ Z. (Skewed boxes won't satisfy `computeNmatrix`'s ±X/±Y normal assertion.)
- **Ghost nodes outside the domain.** For a LEFT (−X) face the ghost column is at
  `x_face − lx`; RIGHT at `x_face + lx`; FRONT/BACK along ∓Y; BOTTOM at `z_face −
  lz`. The element infers FF-vs-SS purely from coordinates + `btype`, so the ghost
  nodes simply have to be on the *outward* side.
- **One ghost layer, conforming.** Don't double-wrap; don't leave the ghost outer
  face tied to anything but the absorbing elements.

### 9.3 `btype` decision table

For each ghost brick, OR together the faces of the **soil domain** it abuts:

| Faces the brick touches | `btype` (2D) | `btype` (3D) | element count |
|---|---|---|---|
| left side | `"L"` | `"L"` | 1 |
| right side | `"R"` | `"R"` | 1 |
| front / back (3D) | — | `"F"` / `"K"` | 1 |
| bottom | `"B"` | `"B"` | 1 |
| bottom + a vertical (edge) | `"BL"`/`"BR"` | `"BL"`/`"BR"`/`"BF"`/`"BK"` | **1** |
| two verticals (vertical edge, 3D) | — | `"LF"`/`"LK"`/`"RF"`/`"RK"` | **1** |
| bottom + two verticals (corner, 3D) | — | `"BLF"`/`"BLK"`/`"BRF"`/`"BRK"` | **1** |

The 2D corners are the bottom-edge cases (`"BL"`, `"BR"`). 3D has the full
side/edge/corner taxonomy (the element's 17 combinations). **Always one element per
brick with the combined string** — the C++ does the tributary weighting.

### 9.4 The emission + driving sequence

apeGmsh's deck/bridge should emit, in order:

1. ghost nodes (`node`), with the soil boundary nodes already present;
2. `element ASDAbsorbingBoundary2D/3D ...` per ghost brick (stage 0 implicitly);
3. `region <tag> -ele <absorb_tags> -rayleigh <aM> 0 <bK> 0` (free-field damping);
4. the base **velocity** `timeSeries` + `-fx/-fy/-fz` on the **bottom** elements
   (or, if using DRM, the DRM ring instead — not both);
5. the **staging hook**: after the gravity `Static` solve + `loadConst`, emit
   `setParameter -val 1 -ele <absorb_tags> stage`, then the `Transient` block.

A natural apeGmsh API shape (sibling of `g.reinforce`, `g.constraints`):

```python
g.absorbing_boundary(
    faces      = <physical group(s) / face set>,   # which domain faces to wrap
    soil       = <material or (G|E, v, rho)>,       # props for the ghost layer
    thickness  = <t>,                               # 2D only
    rayleigh   = (alphaM, betaK),                   # REQUIRED — free-field damping
    base_input = {"x": ts_vx, "z": ts_vz},          # velocity TimeSeries on bottom (optional)
    # generator: extrude conforming ghost layer (match soil element size),
    #            one element/brick, btype from face membership, corners combined.
)
# ... and an analysis-side hook so the driver emits, between gravity and transient:
g.opensees.stage_absorbing()   # -> setParameter -val 1 -ele <tags> stage
```

### 9.5 What to read back

These elements expose only scalar element responses (no stress/strain field):
`eleResponse $tag stage|G|v|rho|E`
([`setResponse`, line 861](SRC/element/absorbentBoundaries/ASDAbsorbingBoundary2D.cpp)).
The physically interesting outputs (boundary reactions, absorbed energy) are read
from the **adjacent soil nodes/elements**, not the absorbing element. So apeGmsh's
results layer needs **no new reader** for these — just:

- a way to **tag** the absorbing elements (so they can be excluded from
  stress/contour plots, where they'd be meaningless), and
- optionally a `stage` probe for diagnostics (confirm the 0→1 switch fired).

### 9.6 Layered & nonlinear soil — what the generator must encode

This is the part most likely to be modelled wrong, so the generator should own it.

**The boundary is LINEAR ELASTIC.** The free-field column and the dashpot
impedances are built from the constant `G, v, rho` you pass — there is **no
constitutive model, no state, no yielding** in the absorbing element (the free-field
stiffness is a constant elastic `C(λ,μ)` re-evaluated each step; the dashpots use
`ρV_p`, `ρV_s` from those constant moduli). The embedded free-field column **never
softens**. Consequences the generator must design around:

1. **Layering → per-element properties.** Each absorbing element carries its own
   `G, v, rho`. Along a side boundary the *stack* of elements reproduces a **layered
   1-D site response** when each element is given the properties of the **layer it
   sits in**. So the generator must resolve, per ghost cell, the local soil
   layer/material and emit that cell's `G, v, rho` — not one global value. (It already
   knows the cell's location and the soil region it abuts.)

2. **Use the SMALL-STRAIN modulus `G₀`,** not a degraded/secant `G`. `V_s = √(G₀/ρ)`.
   When the soil region is given as a nonlinear material, the generator should pull
   the material's **initial** shear modulus (or accept an explicit `G0` override),
   never a back-figured secant value.

3. **Pressure-dependent soil → depth-graded `G₀`.** For `PDMY02`/`PM4Sand`/etc. the
   small-strain stiffness rises with confining stress, hence with **depth**, even
   within one material. The generator must give each side-boundary element a
   **depth-appropriate `G₀(z)`** (treat the depth-graded `G₀` like a finely-layered
   profile). A single bulk `G` down a tall side boundary is wrong. Provide a hook:
   `soil=` may be a callable `G0(x,y,z)` / a per-layer table, not just a scalar.

4. **Keep nonlinearity in the INTERIOR; the boundary far-field must stay ~linear.**
   Because the boundary free-field cannot yield, the truncation must be placed **far
   enough out that the boundary soil stays low-strain**. If the far field is straining
   hard, the (linear) ASD free-field **diverges** from the interior nonlinear soil →
   the boundary stops matching the free-field → spurious reflections. This is a
   genuine limitation vs FLAC/PLAXIS (whose free-field columns run the *same*
   nonlinear law). The generator can't fix it, but it should **document the rule** and
   ideally **warn** if a structure/nonlinear region sits within ~1–2 wavelengths of
   the boundary.

5. **Gravity (stage 0) before the switch — non-negotiable for pressure-dependent
   soil.** The staged workflow (§7.1) is what lets pressure-dependent interior soil
   build its confining stress: gravity in stage 0 (the ASD penalty constraints give
   the standard roller-sides / fixed-base `K₀` condition), then `setParameter
   stage=1`. The generator's staging hook must emit the gravity solve + `loadConst`
   **before** the stage flip. (The ASD free-field stiffness stays frozen at the passed
   `G₀` — fine *because* the boundary soil is meant to stay near small-strain.)

6. **Rayleigh = small-strain damping of the boundary soil** (§7.3), tuned to the
   *small-strain* (low-strain) damping, since the boundary doesn't develop hysteretic
   damping the way the nonlinear interior does.

**API implication.** `g.absorbing_boundary(soil=…)` should accept (a) a scalar
`(G0|E, v, rho)`, (b) a per-layer table keyed on the soil layers, or (c) a callable
`f(x,y,z) → (G0, v, rho)` for pressure-dependent profiles — and resolve it
**per ghost cell**. Defaulting to the adjacent soil region's *initial* modulus is the
right behaviour; a degraded secant is a bug.

### 9.7 Validation hook (recommended)

The cheapest correctness check apeGmsh can automate: build a **single soil column**
with ASD lateral boundaries + a compliant base, drive it with a known velocity
pulse, and confirm the interior reproduces the analytical 1-D free-field transfer
function. This catches the most common user errors at once — **missing Rayleigh on
the free-field** (§7.3), **wrong input quantity/scale** (§7.2), and **wrong/secant
modulus** (§9.6). Extend it to a 2-layer column to also exercise the per-layer
property resolution. See §10.

---

## 10. Verification status & proposed battery

Every structural fact in this guide (class tags, DOF counts, the four mechanisms,
the stage machine, the static condensation, the penalty scaling, the input
convention) is **confirmed by direct citation** to this worktree's `SRC/`.
**Absorption performance is not runtime-verified here** — there are zero bundled
end-to-end examples in the tree.

The natural empirical pin (same treatment the spring guide got,
[[project_zerolength_spring_guide]]) is a 3-test battery against the CPython-3.12
`opensees.pyd` ([[project_opensees_test_env]]):

1. **1-D column self-consistency** — ASD lateral boundaries on a single soil column
   under a vertical S-wave; the interior must reproduce the analytical free-field
   transfer function (verifies mechanisms 1–3 + the Rayleigh-on-FF requirement).
2. **Compliant-base round-trip** — feed a Ricker *velocity* pulse via `-fx`; confirm
   the surface motion = 2× input minus the absorbed reflection (verifies the
   factor-of-2 and the velocity-input convention).
3. **Half-space radiation** — a point load in a 2-D box; measure reflected energy at
   an ASD boundary vs. a bare `LysmerTriangle` vs. a fixed boundary (quantifies the
   free-field column's benefit at oblique incidence).

---

## 11. Related docs

- [[absorbing_boundaries_and_pml_guide]] — the broader survey (PML family, DRM
  input, the explicit-vs-implicit compatibility matrix). This guide is its §3
  deep-dive.
- [[ladruno_apegmsh_contract]] — the apeGmsh-facing feature reference; this
  boundary's generator row points here.
- [[project_robust_central_difference]] / [[project_energy_balance_feature]] — the
  fork's explicit-dynamics work that consumes the dashpot family.
- [[project_zerolength_spring_guide]] — the verified-claim-battery pattern proposed
  in §10; also documents the `ZeroLength`+`Viscous` hand-rolled dashpot alternative.
- SSI/DRM modeling theory: the `quake-research` skill's `references/ssi-drm.md`.

## 12. Maintenance log

- 2026-06-07 — Created. Source-level deep dive on `LysmerTriangle` +
  `ASDAbsorbingBoundary2D/3D`, reading all four `SRC/` files end-to-end. Establishes
  that the ASD elements are a FLAC/PLAXIS-style **free-field boundary + compliant
  base** (free-field column + FF→soil traction + relative-velocity Lysmer dashpots +
  2ρV·v base input), not bare dashpots — correcting/deepening the §3 description in
  [[absorbing_boundaries_and_pml_guide]]. Added the §9 apeGmsh ghost-layer generator
  contract. Absorption performance still unverified — §10 battery proposed.
- 2026-06-07 — Folded a model-building Q&A into §9: added §9.0 (mental model — it's
  an *element* not a material-on-a-brick, no tie/host, one element per ghost *cell*,
  17 `btype` categories, top is free) and §9.6 (**layered & nonlinear soil** — the
  boundary is linear-elastic, so per-element/depth-graded small-strain `G₀`, keep
  nonlinearity in the interior, gravity-before-switch for pressure-dependent soil,
  `soil=` accepts scalar/per-layer/callable). Added gotchas #10 (linear boundary) and
  #11 (element-not-material / no-tie). Same points added to the
  [[ladruno_apegmsh_contract]] generator note.
