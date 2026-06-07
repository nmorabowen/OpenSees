# Absorbing Boundaries & PML in OpenSees — Engineering Reference

A source-accurate survey of every **wave-absorbing boundary** facility shipped in
this OpenSees tree: the **Lysmer dashpot** elements (`LysmerTriangle`,
`ASDAbsorbingBoundary2D/3D`), the **Perfectly Matched Layer** element family
(`PML2D`, `PML3D`, `PML2DVISCOUS`, `PML3DVISCOUS`, and the `PML2D_3/5/12`
research variants), and the **Domain Reduction Method (DRM)** input machinery
that drives waves *into* a truncated domain bounded by these elements.

All file:line citations are against this worktree's `SRC/`. This is a **read /
document** pass — no source was modified.

> **Scope note.** These are all *upstream* OpenSees facilities (not Ladruno
> fork-authored). This guide documents what exists and how to drive it; it is a
> map for anyone building soil–structure-interaction (SSI) / wave-propagation
> models on this fork. Related fork reference: [[zerolength_and_link_springs_guide]]
> (use case #3 there — dashpot absorbing boundaries via `ZeroLength`+`Viscous`).

---

## 1. The three families — TL;DR

| Family | Elements | Dim | Physics | Truncation quality | Stability w/ explicit | Best use |
|---|---|---|---|---|---|---|
| **Lysmer dashpots** | `LysmerTriangle` (tri, 3-node) | 2D/3D surface | Viscous tractions `ρ·Vp·v_n` (P) + `ρ·Vs·v_t` (S) on a boundary facet | First-order; good at normal incidence, leaks at grazing/low-freq | Yes (pure damping, no stiffness in damping mode) | Simple SSI box, quick truncation |
| **ASD absorbing bnd** | `ASDAbsorbingBoundary2D` (quad), `ASDAbsorbingBoundary3D` (brick) | 2D/3D layer | Lysmer–Kuhlemeyer dashpots **+ free-field spring/stiffness + enforced free-field**, staged | First-order + free-field coupling; far better than bare dashpots | Yes (designed for it) | Production SSI, staged gravity→dynamic |
| **PML** | `PML2D`, `PML3D`, `PML2DVISCOUS`, `PML3DVISCOUS` | 2D/3D volume | Complex-coordinate-stretched continuum absorbing layer (Fortran kernels) | Near-reflectionless across freq & incidence (with enough layers) | **No** — implicit, needs effective-stiffness Newmark coupling | High-fidelity SSI, broadband |
| **DRM (input)** | `H5DRMLoadPattern`, `DRMLoadPattern`, `Plane/ParallelDRMInputHandler` | 2D/3D | Effective forces on a 1-element-thick coupling layer (Bielak) | n/a — *injects* the incident field; pairs with one of the above | Either | Seismic input into a truncated domain |

**Mental model.** DRM puts the *incident* wavefield in; the absorbing boundary
(Lysmer / ASD / PML) lets the *scattered* wavefield out. A complete SSI model is
usually **DRM layer + absorbing layer just outside it**.

---

## 2. Class tags & file map

| Element / pattern | Class tag | Header / source | Author(s) |
|---|---|---|---|
| `LysmerTriangle` | `ELE_TAG_LysmerTriangle = 185` | [SRC/element/absorbentBoundaries/LysmerTriangle.cpp](SRC/element/absorbentBoundaries/LysmerTriangle.cpp) | J. A. Abell (UANDES) |
| `ASDAbsorbingBoundary2D` | `ELE_TAG_ASDAbsorbingBoundary2D = 219` | [SRC/element/absorbentBoundaries/ASDAbsorbingBoundary2D.cpp](SRC/element/absorbentBoundaries/ASDAbsorbingBoundary2D.cpp) | M. Petracca, G. Camata (ASDEA) |
| `ASDAbsorbingBoundary3D` | `ELE_TAG_ASDAbsorbingBoundary3D = 220` | [SRC/element/absorbentBoundaries/ASDAbsorbingBoundary3D.cpp](SRC/element/absorbentBoundaries/ASDAbsorbingBoundary3D.cpp) | M. Petracca, G. Camata (ASDEA) |
| `PML2D` | `ELE_TAG_PML2D = 202` | [SRC/element/PML/PML2D.cpp](SRC/element/PML/PML2D.cpp) + `pml_2d.f` | L. Chen, W. Zhang, E. Taciroglu, P. Arduino |
| `PML3D` | `ELE_TAG_PML3D = 201` | [SRC/element/PML/PML3D.cpp](SRC/element/PML/PML3D.cpp) + `pml_3d.f` | W. Zhang, E. Taciroglu, A. Pakzad, P. Arduino |
| `PML2DVISCOUS` | `ELE_TAG_PML2DVISCOUS = 262` | [SRC/element/PML/PML2DVISCOUS.cpp](SRC/element/PML/PML2DVISCOUS.cpp) | A. Pakzad, P. Arduino |
| `PML3DVISCOUS` | `ELE_TAG_PML3DVISCOUS = 271` | [SRC/element/PML/PML3DVISCOUS.cpp](SRC/element/PML/PML3DVISCOUS.cpp) | A. Pakzad, P. Arduino (kernel: W. Zhang, UCLA) |
| `PML2D_3` / `_5` / `_12` | `259` / `260` / `261` | `SRC/element/PML/PML2D_{3,5,12}.cpp` | A. Torino (PhD thesis variants) |
| `LysmerVelocityLoader` | (load class) | [SRC/domain/load/LysmerVelocityLoader.cpp](SRC/domain/load/LysmerVelocityLoader.cpp) | — |
| `H5DRMLoadPattern` | (pattern) | [SRC/domain/pattern/drm/H5DRMLoadPattern.cpp](SRC/domain/pattern/drm/H5DRMLoadPattern.cpp) | J. Abell |
| `DRMLoadPattern` (legacy) | (pattern) | [SRC/domain/pattern/drm/DRMLoadPattern.cpp](SRC/domain/pattern/drm/DRMLoadPattern.cpp) | — |

Build: [SRC/element/absorbentBoundaries/CMakeLists.txt](SRC/element/absorbentBoundaries/CMakeLists.txt)
and [SRC/element/PML/CMakeLists.txt](SRC/element/PML/CMakeLists.txt) (the latter
also compiles the two Fortran kernels via `OPS_ElementFortran`).

Command registration: modern path in
[SRC/interpreter/OpenSeesElementCommands.cpp:831](SRC/interpreter/OpenSeesElementCommands.cpp)
(`functionMap` entries); Tcl path in
[SRC/element/TclElementCommands.cpp](SRC/element/TclElementCommands.cpp);
broker (parallel/clone) in
[SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp](SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp).

> **`element('PML', ...)` aliasing.** The generic `"PML"` command name does **not**
> map to a fixed class. [OpenSeesElementCommands.cpp:606](SRC/interpreter/OpenSeesElementCommands.cpp)
> dispatches on `OPS_GetNDM()`: `ndm==2 → OPS_PML2D()`, otherwise →
> **`OPS_PML3DVISCOUS()`** (not plain `PML3D`). So `element('PML', …)` in 3D
> gives you the *viscous* PML with the mesh-type API (§5.3), while the legacy
> 8-node `PML3D` is only reachable through internal paths. Be explicit: use the
> class names `PML2D` / `PML3DVISCOUS` directly.

---

## 3. Lysmer dashpot family

> **Deep dive lives in [[lysmer_asd_absorbing_boundaries_guide]].** This section is
> the survey-level summary. The formulation-level study — the four mechanisms
> (free-field column + FF→soil traction + **relative-velocity** dashpots + compliant
> base), the stage machine, the 3D static condensation, and the **apeGmsh ghost-layer
> generator contract** — is in that companion guide. Note it **corrects** an
> imprecision below: the ASD elements are a FLAC/PLAXIS-style *free-field boundary*,
> not "Lysmer + free-field spring"; the dashpots act on `V_ff − V_ss` (so only the
> scattered field is absorbed), the base input is a **velocity** history at `2ρV·v`,
> and the free-field column needs its **own Rayleigh damping**.

### 3.1 `LysmerTriangle` — the bare dashpot facet

3-node triangular facet, 3 DOF/node (9 DOF total), placed on the truncation
surface. Applies Lysmer–Kuhlemeyer viscous tractions proportional to the
boundary-normal and tangential velocities.

**Command:**
```tcl
element LysmerTriangle $eleTag $iNode $jNode $kNode  $rho $Vp $Vs  <$length> <$stage>
```
- `$rho` — mass density of the adjacent soil
- `$Vp`, `$Vs` — P- and S-wave velocities of the adjacent soil
- `$length` *(optional)* — explicit element length scale; if omitted the
  element computes its facet area for the dashpot coefficient
- `$stage` *(optional, default 0)* — operating mode (see below)

Parser: [LysmerTriangle.cpp:67](SRC/element/absorbentBoundaries/LysmerTriangle.cpp)
(`Want: element LysmerTriangle eleTag? iNode? jNode? kNode? rho Vp Vs? <length> <stage>`).

**Stage semantics** (the `stage` member, [LysmerTriangle.cpp:312](SRC/element/absorbentBoundaries/LysmerTriangle.cpp)):

| stage | Behaviour |
|---|---|
| **0** (default) | **Damping mode** — contributes only a `getDamp()` dashpot matrix; no stiffness. Use during dynamics. |
| 1 | **Stiffness mode** — contributes a spring `getTangentStiff()` (uses `element_length`/`L`), no damping. Used to hold the boundary during static/gravity. |
| 2 | Both stiffness **and** damping active. |
| 3 | Reaction-passthrough — whatever is in the internal spring-force vector is added directly as a reaction. |

This is the explicit knob that lets you run **gravity with the boundary held
(stage 1/2), then switch to pure absorbing (stage 0)** — though for
`LysmerTriangle` you set the stage at construction; the ASD elements (§4) do the
switch live via `setParameter`.

**Velocity input** — the companion `LysmerVelocityLoader`
([SRC/domain/load/LysmerVelocityLoader.cpp](SRC/domain/load/LysmerVelocityLoader.cpp),
ctor `(tag, eleTag, dir)`) is an `ElementalLoad` that feeds a prescribed boundary
velocity into a Lysmer element in direction `dir` (0/1/2) — i.e. a way to do
"dashpot + incident velocity" boundary input without the full DRM machinery.

### 3.2 `ASDAbsorbingBoundary2D` / `3D` — Lysmer + free-field, staged

ASDEA's production-grade absorbing layer. Beyond bare dashpots it adds a
**free-field column** (spring stiffness + enforced free-field motion) so that the
boundary both absorbs outgoing waves *and* reproduces the correct free-field
soil response — a big accuracy improvement over `LysmerTriangle` at low frequency
and grazing incidence.

**2D — 4-node quad** ([ASDAbsorbingBoundary2D.cpp:142](SRC/element/absorbentBoundaries/ASDAbsorbingBoundary2D.cpp)):
```tcl
element ASDAbsorbingBoundary2D $tag $n1 $n2 $n3 $n4  $G $v $rho $thickness $btype  <-fx $tsxTag> <-fy $tsyTag>
```
**3D — 8-node brick** ([ASDAbsorbingBoundary3D.cpp:482](SRC/element/absorbentBoundaries/ASDAbsorbingBoundary3D.cpp)):
```tcl
element ASDAbsorbingBoundary3D $tag $n1 ... $n8  $G $v $rho $btype  <-fx $tsx> <-fy $tsy> <-fz $tsz>
```
- `$G` — shear modulus, `$v` — Poisson's ratio, `$rho` — density of the soil
- `$thickness` *(2D only)* — out-of-plane thickness
- `$btype` — **boundary-position string**, telling the element which face of the
  domain it sits on so it can orient the dashpots/free-field correctly:
  - 2D: `"B"` (bottom), `"L"` (left), `"R"` (right) — combinable (e.g. corner)
  - 3D: `"B"`, `"L"`, `"R"`, `"F"` (front), `"K"` (back)
  - internally bit-flagged (`BND_BOTTOM`, `BND_LEFT`, … [ASDAbsorbingBoundary2D.cpp:50](SRC/element/absorbentBoundaries/ASDAbsorbingBoundary2D.cpp))
- `-fx/-fy/-fz $tsTag` *(optional)* — `TimeSeries` tags providing a prescribed
  boundary **velocity/force action** per direction (incident-field input).

**Staged analysis (the key usage pattern).** The element has two stages
([ASDAbsorbingBoundary2D.h:40](SRC/element/absorbentBoundaries/ASDAbsorbingBoundary2D.h)):

```
Stage_StaticConstraint = 0   // gravity / static: boundary acts as a constraint (stiffness)
Stage_Absorbing        = 1   // dynamics: dashpots + free-field active
```

You switch **live**, after gravity, via a parameter named `"stage"`
([ASDAbsorbingBoundary2D.cpp:789](SRC/element/absorbentBoundaries/ASDAbsorbingBoundary2D.cpp)):

```tcl
# build model with ASDAbsorbingBoundary* elements (start in stage 0)
# ... apply gravity, run static analysis ...
setParameter -val 1 -ele $tag1 $tag2 ... stage   ;# flip every absorbing element to Stage_Absorbing
# ... now run the transient (dynamic) analysis ...
```

The transition is **one-way**: 0→1 only; trying to go back, or jumping to any
value other than 1, errors out
([updateParameter, ASDAbsorbingBoundary2D.cpp:816](SRC/element/absorbentBoundaries/ASDAbsorbingBoundary2D.cpp)).
You can also query the current stage through `getResponse("stage")`
([ASDAbsorbingBoundary2D.cpp:910](SRC/element/absorbentBoundaries/ASDAbsorbingBoundary2D.cpp)).

> **Node-DOF requirements.** 2D nodes need `ndm==2` and at least 2 DOF; 3D nodes
> need `ndm==3` and ≥3 DOF (checked in `setDomain`,
> [ASDAbsorbingBoundary2D.cpp:343](SRC/element/absorbentBoundaries/ASDAbsorbingBoundary2D.cpp),
> [ASDAbsorbingBoundary3D.cpp:717](SRC/element/absorbentBoundaries/ASDAbsorbingBoundary3D.cpp)).
> The 3D element also rejects a singular Jacobian (excessively distorted brick,
> [ASDAbsorbingBoundary3D.cpp:397](SRC/element/absorbentBoundaries/ASDAbsorbingBoundary3D.cpp)).

**Choosing between the two Lysmer-family elements.** Prefer
`ASDAbsorbingBoundary2D/3D` for any real SSI study — the free-field coupling and
the clean staged gravity→dynamic switch make it the modern default.
`LysmerTriangle` is lighter and triangular (handy on unstructured surface meshes
/ as a quick dashpot skin), but it is a bare first-order dashpot with no
free-field correction.

---

## 4. Perfectly Matched Layer (PML) family

A PML is a finite layer of specially-formulated continuum elements wrapping the
truncation surface that absorbs outgoing waves **near-reflectionlessly across
frequency and angle of incidence** — the gold standard for truncating
wave-propagation domains. In OpenSees the heavy math lives in **Fortran kernels**
(`pml_2d.f`, `pml_3d.f`) wrapped by thin C++ elements that assemble the
mass/damping/stiffness and the PML's extra internal (split-field / stretching)
DOFs.

### 4.1 DOF & wrapper structure

| Element | Nodes | DOF/element | Props | Kernel |
|---|---|---|---|---|
| `PML2D` | 4 | 20 (`PML2D_NUM_DOF`) | 11 (`PML2D_NUM_PROPS`) | `pml_2d.f` |
| `PML3D` | 8 | 72 (`PML3D_NUM_DOF`) | 12 (`PML3D_NUM_PROPS`) | `pml_3d.f` |
| `PML2DVISCOUS` | 4 | 20 | 11 | `pml_2d.f` (+ viscous) |
| `PML3DVISCOUS` | 8 | 72 | 12 | `pml_3d.f` (+ viscous) |

The high DOF counts (20 for a 4-node 2D element, 72 for an 8-node brick) reflect
the PML's **auxiliary stretched-coordinate DOFs** carried per node on top of the
displacement DOFs. From [PML2D.h:43](SRC/element/PML/PML2D.h),
[PML3D.h:43](SRC/element/PML/PML3D.h).

### 4.2 `PML2D` / `PML3D` — the classic elements

`PML2D` parser ([PML2D.cpp:53](SRC/element/PML/PML2D.cpp)):
```tcl
element PML2D $eleTag  $n1 $n2 $n3 $n4   <11 material/PML properties...>
```
It reads **5 ints** (tag + 4 node tags) then **11 doubles** of properties (the
last two are auto-zeroed if you stop short). The 11 props are the soil
constitutive + PML-stretching parameters consumed by `pml_2d.f`.

`PML3D` parser ([PML3D.cpp:57](SRC/element/PML/PML3D.cpp)):
```tcl
element PML3D $eleTag  $n1 ... $n8   $gamma $beta $eta   <12 material/PML properties...>
```
Note `PML3D` additionally pulls **3 Newmark parameters** (`gamma`, `beta`, plus a
time-step term) up front — they are *static* members of the class
([PML3D.cpp:101](SRC/element/PML/PML3D.cpp)) because the PML's effective-stiffness
assembly is integrator-aware. This is the structural tell that **PML is an
implicit-only facility**: it forms an effective stiffness `Keff` tied to the
Newmark coefficients and cannot be driven by a pure explicit central-difference
scheme without that coupling.

### 4.3 `PML3DVISCOUS` — the modern, mesh-type-aware PML

This is what `element('PML', …)` resolves to in 3D (§2). It adds (a) viscous
damping and (b) a convenient **mesh-type geometry API** so you don't hand-compute
the PML reference point and outward normal per element.

Parser ([PML3DVISCOUS.cpp:290](SRC/element/PML/PML3DVISCOUS.cpp)):
```tcl
element PML3DVISCOUS $eleTag  $n1 ... $n8  $matTag $PMLThickness  $meshType <meshType params...>  <options...>
```
- **`$matTag`** — must reference an **`ElasticIsotropicMaterial`** (E, ν, ρ pulled
  from it; other material types are rejected,
  [PML3DVISCOUS.cpp:345](SRC/element/PML/PML3DVISCOUS.cpp)).
- **`$PMLThickness`** — physical thickness of the PML layer.
- **`$meshType`** + params — how the element locates itself relative to the
  absorbing region ([PML3DVISCOUS.cpp:398](SRC/element/PML/PML3DVISCOUS.cpp)):

| `meshType` | Extra params | Status |
|---|---|---|
| `"General"` | `Xref Yref Zref N1 N2 N3` (reference point + outward normal) | implemented |
| `"Box"` | `XC YC ZC L W H` (box center-at-surface + dims; normal auto-computed via `OPS_PML3DVISCOUS_Box`) | implemented |
| `"Sphere"` | `XC YC ZC R` | **NOT implemented yet** (errors out, [PML3DVISCOUS.cpp:462](SRC/element/PML/PML3DVISCOUS.cpp)) |
| `"Cylinder"` | `XC YC ZC R H dir` | **NOT implemented yet** ([PML3DVISCOUS.cpp:487](SRC/element/PML/PML3DVISCOUS.cpp)) |

- **Trailing options** ([PML3DVISCOUS.cpp:509](SRC/element/PML/PML3DVISCOUS.cpp)):
  - `-Newmark $gamma $beta $dt $alpha_f` — integrator coupling (defaults
    `gamma=0.5, beta=0.25`)
  - `-Cp $Cp` — reference phase velocity for the stretching profile
  - `-m $m` — PML polynomial-profile order
  - `-R $R` — target reflection coefficient (theoretical reflectivity the
    stretching is tuned to)
  - `-alphabeta $alpha0 $beta0` — explicitly set the PML stretch coefficients
    instead of deriving them from `Cp/m/R`

`PML2DVISCOUS` is the 2D analogue (same 4-node/20-DOF footprint as `PML2D`, with
viscous damping added).

### 4.4 `PML2D_3` / `PML2D_5` / `PML2D_12` — research variants

Three alternative 4-node 2D PML formulations derived in **Adriano Torino's PhD
thesis** (tags 259/260/261). They share the `PML2D` footprint but differ in the
internal stretched-field formulation. Treat them as experimental / for
reproducing that thesis; the mainstream choice is `PML2D` / `PML2DVISCOUS`.

### 4.5 PML practical notes

- **Implicit only.** PML couples to Newmark coefficients (§4.2). It is **not**
  compatible with the fork's explicit integrators (`CentralDifferenceLadruno`,
  `ExplicitBathe`, …) — for explicit dynamics use the Lysmer / ASD dashpot
  families instead. (Cross-ref the fork's explicit work:
  [[project_robust_central_difference]], [[project_energy_balance_feature]].)
- **Layer thickness vs. accuracy.** Reflectivity is governed by `PMLThickness`,
  the profile order `m`, and target `R`. Too thin a layer or too few elements
  through the thickness leaks low frequencies — the usual PML tradeoff.
- **Cost.** 72 DOF per brick is heavy; PML is the most expensive truncation
  option but also the most accurate. Reserve the layer for the actual boundary
  ring, not the interior.
- **Material restriction.** `PML3DVISCOUS` only accepts `ElasticIsotropic` soil.
  The PML region must be linear elastic (this is standard — PMLs are derived for
  linear media; nonlinearity belongs in the interior, inside the DRM box).

---

## 5. Domain Reduction Method (DRM) — driving the waves in

The absorbing boundary lets *scattered* energy leave; the **DRM** injects the
*incident* seismic field. It does so with **effective forces on a single layer of
elements** straddling the boundary between the "outside" (free-field) and
"inside" (the reduced domain) — Bielak's two-step method. The incident
displacement/acceleration history is read from a file and converted to nodal
effective forces each step.

Files ([SRC/domain/pattern/drm/](SRC/domain/pattern/drm/)):

| Component | Role |
|---|---|
| `H5DRMLoadPattern` | **Modern** DRM — reads the incident field from an **HDF5** `.h5drm` dataset; the production path. |
| `DRMLoadPattern` / `DRMLoadPatternWrapper` | Legacy DRM load pattern (the `PBowlLoading`-era API, [SRC/doc/DRMCommands.tex](SRC/doc/DRMCommands.tex)). |
| `DRMInputHandler` (base) | Abstract incident-field reader. |
| `PlaneDRMInputHandler` | 2D / plane-wave incident field. |
| `ParallelDRMInputHandler` | MPI-partitioned DRM (large 3D runs). |
| `DRMBoundaryLayerDecorator` | Tags / manages the 1-element-thick DRM coupling layer. |

**How it pairs with the boundary.** A canonical SSI model is, from inside out:
1. **Interior** — your structure + nonlinear soil.
2. **DRM layer** — one element ring carrying the `H5DRMLoadPattern` effective forces (the seismic input).
3. **Absorbing layer** — `ASDAbsorbingBoundary*` (or PML, or `LysmerTriangle`) just outside the DRM layer, absorbing the scattered field.

> **Fork history.** This fork already touched DRM once: a `cFactor*DRM_F` scaling
> fix landed on `ladruno` — see [[project_handoff_jaabell_drm_diff]]. That work
> confirms `H5DRMLoadPattern` is the live path here.

---

## 6. Operational scoping — explicit vs implicit (vs the code)

This section is the **real operational map**: what an actual analysis must look
like for each boundary, separated by **implicit** (Newmark / HHT) and **explicit**
(central-difference) transient analysis, with the constraints read straight from
the element source — not the absorption theory.

> **Up front: there are NO bundled end-to-end examples.** A repo-wide sweep
> (`EXAMPLES/`, `Workshops/`, verification scripts, `Ladruno_implementation/`)
> found **zero** complete scripts driving any of these elements through a full
> analysis. Treat every recipe below as **scoped-from-source, not example-backed**
> — these are research-grade elements and the operational details here are
> reconstructed from the C++/Fortran, so verify against a known benchmark before
> production use.

### 6.1 Master compatibility matrix

| Boundary | Implicit (Newmark) | Explicit (central diff) | Hard operational gate (from code) |
|---|---|---|---|
| `LysmerTriangle` (stage 0) | ✅ | ✅ (its natural home) | Pure `getDamp()` only; **no mass, no stiffness** in stage 0 ([LysmerTriangle.cpp:380](SRC/element/absorbentBoundaries/LysmerTriangle.cpp)). Mass must come from the adjacent soil solids. |
| `ASDAbsorbingBoundary2D/3D` | ✅ (primary) | ⚠️ usable | Supplies **its own** `M`,`C`,`K` in stage 1; **staged** — must flip `stage` 0→1 between gravity and dynamics. `M`/`C` are coupling (non-diagonal) matrices. |
| `PML2D` / `PML2DVISCOUS` | ✅ **only** | ❌ | Newmark γ,β **baked in / static members**; forms `Keff=K+(η·dt/β)·G`. Nodes need **ndf = 5**. |
| `PML3D` / `PML3DVISCOUS` | ✅ **only** | ❌ | Same as above; nodes need **ndf = 9**. |
| `H5DRMLoadPattern` (input) | ✅ | ✅ | Integrator-agnostic — it just adds effective nodal forces each step. |

### 6.2 The hard gates, proven from source

**(a) PML requires extra nodal DOFs — this is the #1 build-breaker.**
`PML3D::update()` reads **9 components per node**
([PML3D.cpp:230-238](SRC/element/PML/PML3D.cpp): `uNode(j)` for `j=0..8`), and
`PML3D_NUM_DOF=72 / 8 nodes = 9`. The 6 extra DOFs beyond the 3 displacements are
the PML's **split-field / stretched-coordinate history** carried *on the nodes*.
So the PML region must be built with **`model -ndm 3 -ndf 9`** (2D PML: 20/4 = **5
dof/node → `-ndf 5`**). There is **no runtime ndf check** — too few DOFs reads
out-of-bounds memory and silently corrupts results. At the PML↔soil interface the
shared nodes carry the full ndf=9/5; the regular soil elements only bind their
first 2/3 DOFs.

**(b) PML's Newmark γ,β are baked in and globally shared.**
They are `static` class members set at construction
([PML3D.cpp:101,150-153](SRC/element/PML/PML3D.cpp)) — *all* PML elements share
one γ,β, and **they must equal the integrator's** `Newmark $gamma $beta`. The
effective stiffness `Keff = K + (η·dt/β)·G` is re-formed every
`getTangentStiff()` ([PML3D.cpp:250-257](SRC/element/PML/PML3D.cpp)) and `dt` is
refetched per step (`Domainptr->getDT()`, [PML3D.cpp:222](SRC/element/PML/PML3D.cpp)),
so a **varying dt is tolerated**, but γ/β/integrator mismatch silently gives the
wrong solution. This Newmark coupling is exactly why **PML is implicit-only** —
there is no diagonal-mass / explicit path.

**(c) PML matrices are symmetric → a symmetric solver is fine.**
The code comments `k and g are symmetric matrices`
([PML3D.cpp:252](SRC/element/PML/PML3D.cpp)); nothing forces a nonsymmetric
solver. `system UmfPack` / `BandSPD` / `ProfileSPD` all work — no `FullGeneral`
needed.

**(d) ASD boundary is genuinely staged in the matrices.**
Stage 0 (`Stage_StaticConstraint`) contributes only a **penalty stiffness**
(`addKPenaltyStage0`, [ASDAbsorbingBoundary2D.cpp:444](SRC/element/absorbentBoundaries/ASDAbsorbingBoundary2D.cpp))
— it *holds* the boundary during gravity. Stage 1 (`Stage_Absorbing`) swaps that
for free-field **damping** (`addCff`+`addClk`, line 470) and **mass** (`addMff`,
line 486). It manages its **own** absorbing damping — it does **not** ride on
global Rayleigh. Nodes need ndf=2 (2D) / 3 (3D), standard solid DOFs.

**(e) Lysmer stage 0 is a bare dashpot.**
Stiffness only assembles for `stage∈{1,2}`
([LysmerTriangle.cpp:312](SRC/element/absorbentBoundaries/LysmerTriangle.cpp));
damping `C=diag(ρVs,ρVs,ρVp)` assembles for `stage∈{0,2,3}`
([LysmerTriangle.cpp:380](SRC/element/absorbentBoundaries/LysmerTriangle.cpp)).
No mass is provided in any stage.

### 6.3 Implicit recipe (the supported path for everything)

**PML (3D, `PML3DVISCOUS`) — implicit Newmark:**
```tcl
model BasicBuilder -ndm 3 -ndf 9          ;# <-- 9 DOF for PML-region nodes (3D)
# ... build interior soil (ndf 3 elements) + PML brick ring ...
element PML3DVISCOUS $tag $n1 .. $n8  $matTag $thick  Box  $XC $YC $ZC $L $W $H  -Newmark 0.5 0.25 $dt 0.0
# ... constraints, loads ...
constraints Transformation
numberer RCM
system UmfPack                            ;# symmetric OK (§6.2c)
algorithm Linear                          ;# PML region is linear-elastic
integrator Newmark 0.5 0.25               ;# gamma,beta MUST match the element's
analysis Transient
analyze $nsteps $dt                       ;# keep dt consistent with the element
```

**ASD absorbing boundary — implicit, staged gravity→dynamics:**
```tcl
model BasicBuilder -ndm 3 -ndf 3
# ... soil + ASDAbsorbingBoundary3D elements (created in stage 0) ...
# --- gravity (stage 0: boundary held by penalty) ---
integrator LoadControl 1.0
analysis Static; analyze 1
loadConst -time 0.0
# --- flip every boundary element to absorbing mode ---
setParameter -val 1 -ele $b1 $b2 ... stage      ;# 0 -> 1, one-way
# --- dynamics (stage 1: free-field M/C/K active) ---
integrator Newmark 0.5 0.25
analysis Transient; analyze $nsteps $dt
```

### 6.4 Explicit recipe (Lysmer / ASD only — PML excluded)

For explicit central-difference the boundary must not inject stiffness that
forces a solve, and the global mass should be (near-)diagonal. **PML cannot
participate** (§6.2b). The dashpot families can:

```tcl
model BasicBuilder -ndm 3 -ndf 3
# ... soil solids supply the (lumped) mass; Lysmer facets supply boundary damping ...
element LysmerTriangle $tag $i $j $k  $rho $Vp $Vs        ;# stage 0 = pure dashpot
# (or ASDAbsorbingBoundary3D, already in absorbing stage for the dynamic run)
system Diagonal                          ;# explicit wants a diagonal/lumped system
algorithm Linear
integrator CentralDifference             ;# or the fork's CentralDifferenceLadruno
analysis Transient
analyze $nsteps $dt                      ;# dt < dt_cr (CFL of the soil mesh)
```

> **Explicit caveat (operational, not method):** `LysmerTriangle` stage-0 damping
> is **diagonal** (`diag(ρVs,ρVs,ρVp)`), so it folds cleanly into an explicit
> step. `ASDAbsorbingBoundary*` supplies **coupling (non-diagonal) `M`/`C`** at the
> boundary ring, so those boundary DOFs are not strictly diagonal — workable, but
> it is the reason `LysmerTriangle` is the cleaner *explicit* dashpot and ASD is
> the cleaner *implicit, free-field-correct* boundary. For explicit work prefer
> the fork's `CentralDifferenceLadruno` ([[project_robust_central_difference]]),
> noting that **damping shrinks `dt_cr`** in explicit schemes.

### 6.5 What this means for use-case selection

- **Implicit SSI, broadband fidelity, linear far-field** → PML (mind ndf=9/5,
  matched Newmark, symmetric solver). Highest accuracy, implicit cost.
- **Implicit SSI, gravity-then-quake, nonlinear interior** → ASD absorbing
  boundary with the staged switch. The pragmatic production default.
- **Explicit dynamics (impact, blast, high-rate, fork explicit integrators)** →
  `LysmerTriangle` stage-0 dashpots (mass from soil), or ASD in absorbing stage.
  **PML is off the table.**
- **Any of the above + seismic input** → `H5DRMLoadPattern` ring inside the
  boundary; integrator-agnostic.

---

## 7. Decision guide

**"Which boundary do I use?"**

- **Quick truncation, explicit dynamics, unstructured surface** → `LysmerTriangle`
  (stage 0). Cheapest, first-order.
- **Production SSI, gravity-then-earthquake, want free-field correction** →
  `ASDAbsorbingBoundary2D/3D` with the staged `setParameter ... stage` switch.
  The pragmatic default. Works explicit or implicit.
- **Broadband, high-fidelity, low-frequency-critical, implicit run** → **PML**
  (`PML3DVISCOUS` in 3D via the `Box`/`General` mesh-type API; `PML2DVISCOUS` in
  2D). Most accurate, most expensive, implicit-only, linear-elastic layer.
- **Seismic input into any of the above** → `H5DRMLoadPattern` on a one-element
  DRM ring just inside the absorbing layer.

**"Explicit or implicit?"**

| Boundary | Explicit OK? |
|---|---|
| `LysmerTriangle` (stage 0) | ✅ pure damping matrix |
| `ASDAbsorbingBoundary*` | ✅ designed for it |
| `PML*` | ❌ implicit (Newmark-coupled effective stiffness) |

---

## 8. Gotchas worth pinning

0. **PML nodes need extra DOFs — `-ndf 9` (3D) / `-ndf 5` (2D).** The auxiliary
   split-field DOFs are real nodal DOFs with **no runtime check**; too few = silent
   out-of-bounds. The single most common PML build-breaker (§6.2a).
1. **`element('PML', …)` is NDM-dispatched** and in 3D gives `PML3DVISCOUS`, not
   the classic `PML3D`. Name the class explicitly (§2).
2. **ASD absorbing boundaries are staged and the switch is live, one-way.** Build
   in stage 0, run gravity, then `setParameter -val 1 -ele <tags> stage` before
   the dynamic step. You cannot go back to 0 (§3.2).
3. **`$btype` strings orient the element.** Getting `"L"/"R"/"B"/"F"/"K"` wrong
   (or omitting the corner combination) silently mis-orients the dashpots /
   free-field. Match the string to the physical face.
4. **PML soil must be linear elastic** — `PML3DVISCOUS` rejects anything but
   `ElasticIsotropicMaterial`. Put nonlinearity inside the DRM box, never in the
   PML ring.
5. **PML `Sphere`/`Cylinder` mesh types are stubs** — they parse then error. Only
   `General` and `Box` actually build a PML (§4.3).
6. **PML is implicit.** Don't pair it with the fork's explicit integrators.
7. **DRM + absorbing layer are complementary, not alternatives.** DRM injects;
   the boundary absorbs. A model that has one without the other is usually wrong.

---

## 9. Verification status

Source-level facts in this guide (class tags, DOF counts, parser argument lists,
the `stage` enum and one-way `updateParameter` guard, the `meshType` switch and
its unimplemented branches, the `element('PML')` NDM dispatch) are **confirmed by
direct citation** to this worktree's `SRC/`. Quantitative absorption performance
(reflectivity vs. layer count, low-frequency leakage) is **not** runtime-verified
here — that would need a dedicated wave-propagation test battery (a half-space
under a Ricker pulse, measuring reflected energy at the boundary), which is a
natural follow-up if we want the same empirical-pin treatment the spring guide
([[zerolength_and_link_springs_guide]]) got.
