# ZeroLength & Link "Spring" Elements — Engineering Reference

A deep, source-accurate guide to the OpenSees "spring" element family:
**ZeroLength**, **ZeroLengthSection**, **ZeroLengthND**, the **ZeroLengthContact**
variants, **CoupledZeroLength**, and **TwoNodeLink** (the "two-point spring").

The focus is **stiffness AND damping behavior**, and suitability for four use
cases: (1) contact, (2) phenomenological hinges, (3) absorbing boundaries
(dashpots), (4) penalty/constraint emulation.

All file:line citations are against this worktree's `SRC/element/`.

> **Verification status.** Source-level facts (defaults, `setTrialStrain`
> signatures, assembly) are confirmed by direct citation. The *surprising
> behavioural* claims are additionally **pinned by a runtime battery** —
> `tests/test_spring_damping_claims.py` (8/8 passing, CPython 3.12 + built
> `dist/bin/opensees.pyd`) — plus the pre-existing `tests/test_damping_channels.py`
> for the plain-ZeroLength `-doRayleigh` trap. Empirically verified here:
> ZeroLengthSection `-doRayleigh` default **ON** (§2.2/§3.2); ZeroLength+`Viscous`
> damps with no flag (§4.3); `Viscous` alone gives a singular static tangent
> (§5.2); the **double-negative** — `Viscous` is inert inside
> SectionAggregator→ZeroLengthSection (§2.7/§4.2); TwoNodeLink applies material
> damping with no flag (§2.5); CoupledZeroLength yields on the resultant (§2.4).
> Testing also surfaced the **ndf=3/6 requirement** for ZeroLengthSection (§2.2).

---

## 1. Overview / TL;DR table

| Element | Materials taken | DOFs handled | Carries mass? | Velocity-dependent (dashpot) material works? | Rayleigh damping default | Best use |
|---|---|---|---|---|---|---|
| **ZeroLength** | N × `uniaxialMaterial`, each on a chosen local dir 1–6 | 1/2/3D, 2–12 numDOF | No | **YES** — `setTrialStrain(strain, strainRate)` ([ZeroLength.cpp:819](#)) | **OFF** (`-doRayleigh` default 0) | Hinges, dashpots, penalty springs, gap springs |
| **ZeroLengthSection** | 1 × `section` (fiber/aggregator) | 3 (2D) or 6 (3D) per node | No | **NO** — section deformation only, no rate ([ZeroLengthSection.cpp:289](#)) | **ON** (`-doRayleigh` default 1) | Fiber lumped-plastic hinge with axial-moment coupling |
| **ZeroLengthND** | 1 × `nDMaterial` (order 2/3/5/6) + optional 1 × `uniaxialMaterial` | 3 (2D) or 6 (3D) per node | No | **NO** — `setTrialStrain` rate-less ([ZeroLengthND.cpp:383](#)) | none (always zero) | Coupled multi-axial constitutive hinge (e.g. soil cap) |
| **CoupledZeroLength** | 1 × `uniaxialMaterial` acting on the **resultant** of 2 dirs | 2D/3D | No | **YES** — `setTrialStrain(strain, strainRate)` ([CoupledZeroLength.cpp:435](#)) | **OFF** (`<useRayleigh?>` default 0) | Bi-directional (kinematic) coupled spring/bearing |
| **TwoNodeLink** | N × `uniaxialMaterial`, each on a chosen local dir 1–6 | 1/2/3D, 2–12 numDOF | **Yes** (`-mass`, lumped) | **YES** — `setTrialStrain(ub, ubdot)` ([TwoNodeLink.cpp:663](#)) | **OFF** (`-doRayleigh`) **but material damping always on** | Bearings, dampers, P-Delta links, dashpots with mass |
| **ZeroLengthContact2D/3D** | none — built-in penalty (Kn, Kt, fs) | 2D (2 dof/node) / 3D (3 dof/node) | No | No (penalty contact, no rate) | none | Node-to-node penalty contact + Coulomb friction |
| **ZeroLengthContactASDimplex** | none — Kn, Kt, mu, IMPL-EX | 2D/3D | No | No (rate-independent friction via IMPL-EX) | none | Robust explicit-friendly node-to-node contact |
| **ZeroLengthImpact3D** | none — Hertz-damp impact params | 3D | No | Impact damping is internal | none | Pounding/impact |

**The single most important gotcha:** the `-doRayleigh` default is *not*
consistent across the family. It defaults **OFF** for ZeroLength / CoupledZeroLength
/ TwoNodeLink, but **ON** for ZeroLengthSection. See §3.

### Section builders that the spring elements consume

`ZeroLengthSection` (and force/displacement beam elements) take a
`SectionForceDeformation`, not raw uniaxial materials. The standard way to
assemble a **multi-DOF uncoupled "section" from uniaxial springs** — and the
canonical companion to `ZeroLengthSection` for phenomenological hinges — is
`section Aggregator`. It is a *section*, not an element.

| Section builder | What it is | Coupling | Strain rate to mats? | Best use |
|---|---|---|---|---|
| **SectionAggregator** (`section Aggregator`) | Glues N `uniaxialMaterial`s onto section response codes (P, Mz, My, Vy, Vz, T), optionally on top of a base `-section` | **Uncoupled — block diagonal** ([SectionAggregator.cpp:496-497](#)) | **NO** — `setTrialStrain(e)` only, rate defaults 0 ([SectionAggregator.cpp:447](#)) | Build a full P-M-V-T hinge section for ZeroLengthSection |

See §2.7 for the deep dive and §4.2 for the canonical hinge wiring.

### The constitutive half — and the richer cousins

A spring element is only a frame; the **material** decides whether it behaves as
contact, a dashpot, or a hinge. **§5 catalogs the spring-defining materials**
(gap/contact, viscous, degrading-hinge) with the all-important *rate-dependent?*
column. **§6 covers the bearing/isolator element family** (elastomeric +
friction-pendulum) — purpose-built coupled springs with a true bidirectional
shear yield surface that generic springs cannot reproduce. **§7 covers the joint
element family** (panel-zone spring bundles) — the specialized cousin of a
phenomenological hinge for when the hinge *is* a beam-column joint.

### Related families — present in this tree but out of scope here

These solve the same use cases with dedicated elements rather than generic
springs, and are *not* deep-dived in this guide (pointers only):
- **Absorbing boundaries:** `ASDAbsorbingBoundary2D/3D`, `LysmerTriangle`
  (+`LysmerVelocityLoader`), and the **PML** family (`PML2D/3D`, `PML2D_3/5/12`,
  `PML*VISCOUS`) — the engineered alternative to a ZeroLength+`Viscous` dashpot
  (§4.3).
- **Surface / embedded contact & ties:** `ASDEmbeddedNodeElement`,
  `EmbeddedBeamInterface*` (non-matching-mesh node ties / typed constraints).
- **Typed constraints proper:** `equalDOF`, `rigidLink`, `rigidDiaphragm`
  (MP_Constraints) — the exact-kinematic alternative to penalty springs (§4.4).

---

## 2. Per-element deep dives

### 2.1 ZeroLength

**Source:** `SRC/element/zeroLength/ZeroLength.{h,cpp}`. Parser:
`OPS_ZeroLength()` ([ZeroLength.cpp:61](#)).

#### Command signature

```
element zeroLength $eleTag $iNode $jNode \
    -mat $matTag1 $matTag2 ... \
    -dir $dir1 $dir2 ... \
    <-orient $x1 $x2 $x3 <$yp1 $yp2 $yp3>> \
    <-doRayleigh $rFlag> \
    <-dampMats $dmpTag1 $dmpTag2 ...> \
    <-damp $dampingTag>
```

- **Materials:** any number of `uniaxialMaterial`s. Parsed greedily after `-mat`
  ([ZeroLength.cpp:107-121](#)). Each material maps to exactly one local direction.
- **`-dir` (alias `-dof`):** one direction per material, 1–6 = local
  dx, dy, dz, rx, ry, rz. Internally decremented to 0–5
  ([ZeroLength.cpp:162-164](#)); `checkDirection` clamps out-of-range to 0
  ([ZeroLength.cpp:1891-1898](#)). Special case: a `dir==2` in a 2D model is
  silently remapped to `5` ("For Keri Ryan", [ZeroLength.cpp:297-299](#)) so that
  "2" means rotation `Mz` in 2D-with-3-dof.

#### Direction → DOF / orientation semantics

The element is a **relative-deformation spring**: deformation for material *m*
is `eps_m = (u_J − u_I) projected onto local axis indx`, where `indx = dir % 3`
selects the translation/rotation axis ([ZeroLength.cpp:1927-1930](#)).

`-orient x1 x2 x3 yp1 yp2 yp3` defines the local frame: local **x** = given x,
**z = x × yp**, **y = z × x**, all normalized into `transformation(3,3)`
([ZeroLength.cpp:1857-1884](#)). In 2D only x is needed; y is derived as
`(-x1, x0)` ([ZeroLength.cpp:247-251](#)). **Because the nodes are coincident the
geometry carries no orientation information — the local frame comes *entirely*
from `-orient` (defaulting to global x/y).**

#### Stiffness assembly

`setTran1d()` builds the strain-displacement row vector `t1d(mat, dof)`
([ZeroLength.cpp:1903-1995](#)). For a translational dir it places the relevant
row of `transformation` into the node-J half of the DOFs, then mirrors it with a
negative sign into the node-I half ([ZeroLength.cpp:1991-1992](#)). So each
material's row is `t = [−c | +c]` where `c` is the direction cosine vector — the
classic relative-DOF coupling.

`getTangentStiff()` ([ZeroLength.cpp:837-872](#)) forms
`K = Σ_m t1d[m]^T · E_m · t1d[m]`, where `E_m = theMaterial1d[m]->getTangent()`.
`getInitialStiff()` is identical with `getInitialTangent()`
([ZeroLength.cpp:875-910](#)).

#### Damping behavior (important)

`update()` ([ZeroLength.cpp:792-835](#)) computes **both** strain and strain rate
from nodal disp and **vel**, and calls
`theMaterial1d[mat]->setTrialStrain(strain, strainRate)`
([ZeroLength.cpp:817-819](#)). **This is the key enabler for dashpots:** a
`Viscous`, `ViscousDamper`, or any rate-dependent uniaxial material gets a real
strain rate and produces a velocity-proportional force.

`getDamp()` ([ZeroLength.cpp:913-968](#)) has three branches keyed off
`useRayleighDamping`:
- `== 1`: returns `Element::getDamp()` → Rayleigh `αM + βK` matrix.
- `== 2`: uses a parallel set of damping materials' `getTangent()` (the
  `-dampMats` feature).
- **else (`== 0`, the default):** uses each material's `getDampTangent()`
  ([ZeroLength.cpp:952](#)).

> **Trap:** with the default `-doRayleigh 0`, the element's `getDamp()` returns
> `Σ t^T · getDampTangent() · t`. For most uniaxial materials `getDampTangent()`
> returns **0**, so **stiffness-proportional Rayleigh damping contributes nothing
> unless you pass `-doRayleigh 1`.** A `Viscous` material, by contrast, returns a
> nonzero `getDampTangent()` and works regardless. The default is OFF
> ([ZeroLength.cpp:171](#), `int doRayleighDamping = 0;`).

`getResistingForceIncInertia()` ([ZeroLength.cpp:1051-1077](#)) adds Rayleigh
damping forces only if `useRayleighDamping==1` and an α/β is nonzero. The element
has **no mass** (`getMass()` returns zero, [ZeroLength.cpp:971-977](#);
`addInertiaLoadToUnbalance` is a no-op, [ZeroLength.cpp:994-999](#)).

There is also a new `-damp $dampingTag` path that attaches a `Damping` object and
scales stiffness via `getStiffnessMultiplier()` ([ZeroLength.cpp:855](#),
[ZeroLength.cpp:1026-1047](#)).

#### State / commit

The element holds essentially no history of its own except optional initial
offset vectors `d0`/`v0` captured at `setDomain` (so a prestrained/prevelocity
state is subtracted, [ZeroLength.cpp:692-698](#), [ZeroLength.cpp:805-809](#)).
All path-dependent state lives in the uniaxial materials; `commitState` just
forwards to them ([ZeroLength.cpp:733-753](#)).

---

### 2.2 ZeroLengthSection

**Source:** `ZeroLengthSection.{h,cpp}`. Parser: `OPS_ZeroLengthSection()`
([ZeroLengthSection.cpp:53](#)).

#### Command signature

```
element zeroLengthSection $eleTag $iNode $jNode $secTag \
    <-orient $x1 $x2 $x3 <$yp1 $yp2 $yp3>> \
    <-doRayleigh $rFlag>
```

Takes exactly **one `section`** (`SectionForceDeformation`), e.g. a fiber
section or `section Aggregator` ([ZeroLengthSection.cpp:114](#)).

> **Hard requirement (verified by test):** the model **must use `-ndf 3` (2D) or
> `-ndf 6` (3D)** — `setDomain()` aborts with *"element only works for 3 (2d) or
> 6 (3d) dof per node"* ([ZeroLengthSection.cpp:247](#)) otherwise. Unlike plain
> `ZeroLength` (which works at any ndf), ZeroLengthSection always operates on the
> full section-DOF set, so you cannot use it in a reduced `-ndf 2` model. A
> too-low ndf makes the element **silently absent** (the command errors but
> analysis proceeds with no spring). This bit the verification battery — see
> `tests/test_spring_damping_claims.py`.

#### Section DOF → element DOF

The section's `getType()` code (P, Mz, Vy, My, Vz, T) drives `setTransformation()`
([ZeroLengthSection.cpp:816-917](#)). Each section response component is mapped to
the matching element local DOF via the direction-cosine matrix: `SECTION_RESPONSE_P`
→ axial (local x), `MZ` → bending about local z, `T`/`MY`/`VZ` for 3D, etc.
The resulting `A` matrix (order × numDOF) again places `+` on node J and `−` on
node I ([ZeroLengthSection.cpp:913-915](#)).

`getTangentStiff()` = `A^T · kb · A` with `kb = theSection->getSectionTangent()`
([ZeroLengthSection.cpp:325-341](#)); `getResistingForce()` = `A^T · q`
([ZeroLengthSection.cpp:388-404](#)). The big advantage over multi-material
ZeroLength: the section can couple axial and bending (fiber coupling of P–M),
which independent uniaxial springs cannot.

#### Damping behavior

`getDamp()` ([ZeroLengthSection.cpp:343-351](#)):
- `useRayleighDamping == 1` → `Element::getDamp()` (Rayleigh).
- else → **zero matrix**.

`update()` ([ZeroLengthSection.cpp:282-295](#)) calls only
`setTrialSectionDeformation(*v)` — **no strain rate is ever passed**. There is no
section-rate channel. Therefore:

> **ZeroLengthSection cannot act as a dashpot / velocity-dependent element.** Its
> only damping is Rayleigh, and here `-doRayleigh` defaults to **1 / ON**
> ([ZeroLengthSection.cpp:76](#), `int doRayleighDamping = 1;`) — the opposite of
> plain ZeroLength.

No mass. No element-held history beyond the section.

---

### 2.3 ZeroLengthND

**Source:** `ZeroLengthND.{h,cpp}`. Parser: `OPS_ZeroLengthND()`
([ZeroLengthND.cpp:67](#)).

#### Command signature

```
element zeroLengthND $eleTag $iNode $jNode $ndTag <$uniTag> \
    <-orient $x1 $x2 $x3 $yp1 $yp2 $yp3>
```

- One `nDMaterial` of order 2, 3, 5, or 6 ([ZeroLengthND.cpp:155](#)); copied as
  `PlaneStrain2D` (2D) or `ThreeDimensional` (3D) ([ZeroLengthND.cpp:141-144](#)).
- Optional one extra `uniaxialMaterial` that is wired to local-x (transformation
  row index 2, [ZeroLengthND.cpp:1026-1042](#)) to supply, e.g., an axial
  response the ND model lacks.

#### Stiffness / coupling

`getTangentStiff()` = `A^T · kb · A` with the **full (coupled) ND tangent**
`kb = theNDMaterial->getTangent()` — note the double loop over `(k,l)` includes
off-diagonal coupling ([ZeroLengthND.cpp:397-404](#)), unlike the
diagonal-per-material ZeroLength. The optional uniaxial term is added on row 2
([ZeroLengthND.cpp:406-418](#)).

#### Damping behavior

`getDamp()` returns **zero, unconditionally** ([ZeroLengthND.cpp:471-478](#)).
`update`/stiffness call `theNDMaterial->setTrialStrain(*v)` and
`the1DMaterial->setTrialStrain(e)` — **strain only, no rate**
([ZeroLengthND.cpp:383](#), [ZeroLengthND.cpp:409](#)). There is **no
`-doRayleigh` option at all**, and `getResistingForceIncInertia()` simply returns
`getResistingForce()` ([ZeroLengthND.cpp:553-558](#)).

> **ZeroLengthND has no damping path whatsoever** — not even Rayleigh, not even
> via materials. No dashpot capability. No mass.

---

### 2.4 CoupledZeroLength

**Source:** `CoupledZeroLength.{h,cpp}`. Parser:
([CoupledZeroLength.cpp:62](#)).

```
element CoupledZeroLength $tag $iNode $jNode $dirn1 $dirn2 $matTag <$useRayleigh>
```

A **single** uniaxial material acts on the *resultant deformation* of **two**
local directions (`dirn1`, `dirn2`), so the spring is radial/kinematically coupled
in a plane (think a circular-interaction friction or bearing element).
`update()` passes both strain and strain rate:
`theMaterial->setTrialStrain(strain, strainRate)` ([CoupledZeroLength.cpp:435](#))
→ **dashpot-capable**. `getDamp()` mirrors ZeroLength: Rayleigh if `useRayleigh==1`,
else material `getDampTangent()` ([CoupledZeroLength.cpp:530-539](#)).
`<useRayleigh>` defaults **OFF** ([CoupledZeroLength.cpp:81](#)).

---

### 2.5 TwoNodeLink ("two-point spring")

**Source:** `SRC/element/twoNodeLink/TwoNodeLink.{h,cpp}`. Parser:
`OPS_TwoNodeLink()` ([TwoNodeLink.cpp:57](#)). Author: Andreas Schellenberg.

#### Command signature

```
element twoNodeLink $eleTag $iNode $jNode \
    -mat $matTag1 ... -dir $dir1 ... \
    <-orient <$x1 $x2 $x3> $y1 $y2 $y3> \
    <-pDelta $Mratios...> \
    <-shearDist $sDratios...> \
    <-doRayleigh> \
    <-mass $m>
```

- N `uniaxialMaterial`s, one per local direction 1–6 ([TwoNodeLink.cpp:80-99](#)),
  exactly like ZeroLength.
- **Crucially, the element is NOT zero-length-only:** it reads the actual nodal
  geometry to compute its length `L` and local-x from node coordinates
  ([TwoNodeLink.cpp:124-146](#), [TwoNodeLink.cpp:1357-1365](#)). It works fine
  with coincident nodes (then orientation must come from `-orient`), but the
  finite length enables the moment-arm coupling below.

#### Two-step transformation, shearDist and P-Delta (the real difference)

Response is mapped global → local (`Tgl`, direction cosines,
[TwoNodeLink.cpp:1489-1535](#)) → basic (`Tlb`,
[TwoNodeLink.cpp:1540-1575](#)). The basic-system transform is where TwoNodeLink
earns its keep:

- **`-shearDist $sDy <$sDz>`** (ratios in [0,1], default 0.5): a *shear* spring
  (dir 1 or 2) does not just connect the two translational DOFs — it also reacts
  on the rotational DOFs through a moment arm `±sD·L` / `±(1−sD)·L`
  ([TwoNodeLink.cpp:1556-1568](#)). This places the shear-force resultant at a
  user-chosen point along the link, generating the correct end moments. Plain
  ZeroLength cannot do this (its shear and moment springs are independent).
- **`-pDelta $rMy1 $rMy2 $rMz1 $rMz2`**: adds geometric (P-Δ) stiffness and
  forces. The axial force `N` (the dir-0 basic force) multiplies the transverse
  relative displacement to add `N/L` geometric stiffness terms and distribute
  P-Δ moments per the ratios ([TwoNodeLink.cpp:1671-1768](#) for stiffness,
  [TwoNodeLink.cpp:1578-1668](#) for forces). This is what makes TwoNodeLink
  suitable as a *member* (e.g. a base isolator under gravity) rather than a pure
  joint spring.

`getTangentStiff()` = `Tgl^T (Tlb^T kb Tlb + Kpδ) Tgl`
([TwoNodeLink.cpp:669-696](#)); `getInitialStiff` analogous without P-Δ
([TwoNodeLink.cpp:699-721](#)).

#### Damping behavior

`update()` calls `theMaterials[i]->setTrialStrain(ub(i), ubdot(i))`
([TwoNodeLink.cpp:662-663](#)) — strain **and** rate → **dashpot-capable**.

`getDamp()` ([TwoNodeLink.cpp:724-753](#)) is subtly different from ZeroLength:

```
if (addRayleigh == 1)  theMatrix = Element::getDamp();   // Rayleigh part
cb(i,i) = theMaterials[i]->getDampTangent();             // ALWAYS added
theMatrix += Tgl^T (Tlb^T cb Tlb) Tgl;
```

> **Key contrast:** TwoNodeLink **always** adds the materials' `getDampTangent()`
> to the damping matrix, *regardless of* `-doRayleigh`. `-doRayleigh` only
> controls the *additional* Rayleigh `αM + βK` contribution
> ([TwoNodeLink.cpp:731-734](#), default `doRayleigh = 0`,
> [TwoNodeLink.cpp:149](#)). So a `Viscous` material in a TwoNodeLink damps
> whether or not you set `-doRayleigh`.

#### Mass

TwoNodeLink **does** carry mass via `-mass $m` (default 0). It is lumped 50/50 on
the two nodes, **translational DOFs only** ([TwoNodeLink.cpp:756-772](#)), and is
applied in `addInertiaLoadToUnbalance` / `getResistingForceIncInertia`
([TwoNodeLink.cpp:791-818](#), [TwoNodeLink.cpp:845-873](#)).

---

### 2.6 ZeroLengthContact family (summary)

These are **self-contained penalty contact** elements — they take stiffness/
friction *scalars*, not materials, and implement the gap/normal/Coulomb logic
internally. They are the OpenSees answer to unilateral contact between two
coincident (or nearly coincident) nodes.

| Element | Signature | Notes |
|---|---|---|
| **ZeroLengthContact2D** | `element zeroLengthContact2D $tag $iN $jN $Kn $Kt $fs -normal $Nx $Ny` ([ZeroLengthContact2D.cpp:61](#)) | Penalty Kn/Kt + Coulomb `fs`, augmented-Lagrange fields present; 2 dof/node. |
| **ZeroLengthContact3D** | `element zeroLengthContact3D $tag $iN $jN $Kn $Kt $fs $c $dir <$origX $origY>` ([ZeroLengthContact3D.cpp:58](#)) | Adds cohesion `c`; `dir` selects normal axis (0 = circular/auto about an origin). 3 dof/node. |
| **ZeroLengthContactASDimplex** | `element zeroLengthContactASDimplex $tag $iN $jN $Kn $Kt $mu <-orient $x1 $x2 $x3> <-intType $t>` ([ZeroLengthContactASDimplex.cpp:114](#)) | IMPL-EX integration of friction → robust, explicit-friendly; has its own `getDamp()`. The element this fork's Ladruno-contact work builds on. |
| **ZeroLengthContactNTS2D / ZeroLengthInterface2D** | node-to-segment 2D contact / interface | Multi-node, surface-style; secondary. |
| **ZeroLengthImpact3D** | Hertz-contact + impact damping pounding element | Internal impact damping; for pounding studies. |
| **ZeroLengthRocking** | rocking interface | Specialized. |

All of these compute a **contact-state-dependent** tangent (open gap → zero
normal stiffness; closed → Kn; stick/slip → Kt or frictional). None take a
velocity-dependent material; damping (where present, e.g. ASDimplex, Impact3D) is
internal to the element. They are *not* general dashpots.

---

### 2.7 SectionAggregator (the section builder ZeroLengthSection consumes)

**Source:** `SRC/material/section/SectionAggregator.{h,cpp}`. Parser:
`OPS_SectionAggregator()` ([SectionAggregator.cpp:50](#)). Class tag
`SEC_TAG_Aggregator`. Author: MHS (Michael Scott), Jun 2000.

This is **not an element** — it is a `SectionForceDeformation`. It is the standard
way to assemble a multi-DOF, *uncoupled* "section" out of individual
`uniaxialMaterial`s, optionally decorating an existing base section. The original
design intent (header comment, [SectionAggregator.h:32-35](#)) is to "decorate an
MP section (couple bending and axial) with an uncoupled shear relation" — i.e. add
Vy/Vz/T springs to a fiber P-Mz-My section. `ZeroLengthSection` then consumes the
result.

#### Command signature

```
section Aggregator $secTag $matTag1 $code1 $matTag2 $code2 ... <-section $baseSecTag>
```

- **`$matTagI $codeI` pairs:** each `uniaxialMaterial` is glued to one section
  response quantity named by `$codeI` ([SectionAggregator.cpp:98-139](#)). Codes:

  | code string | meaning | enum ([SectionAggregator.cpp:121-132](#)) |
  |---|---|---|
  | `P`  | axial force / axial strain | `SECTION_RESPONSE_P` |
  | `Mz` | bending about local z / curvature z | `SECTION_RESPONSE_MZ` |
  | `My` | bending about local y / curvature y | `SECTION_RESPONSE_MY` |
  | `Vy` | shear along local y / shear strain y | `SECTION_RESPONSE_VY` |
  | `Vz` | shear along local z / shear strain z | `SECTION_RESPONSE_VZ` |
  | `T`  | torsion / twist rate | `SECTION_RESPONSE_T` |

  An unrecognized code aborts parsing with `WARNING invalid code`
  ([SectionAggregator.cpp:133-137](#)).
- **`-section $baseSecTag`** (optional): an existing `SectionForceDeformation`
  (typically a fiber section) whose responses are kept and **prepended**; the
  aggregated uniaxial codes are appended after it ([SectionAggregator.cpp:70-90](#),
  [SectionAggregator.cpp:147-151](#)). The base section is **deep-copied**
  (`theSec.getCopy()`, [SectionAggregator.cpp:228](#)) as are all additions
  (`getCopy()`, [SectionAggregator.cpp:253](#)).
- Related: `section Uniaxial $tag $matTag $code` (`OPS_UniaxialSection()`,
  [SectionAggregator.cpp:156](#)) is the degenerate 1-material case — it builds a
  `SectionAggregator` with no base section and a single addition.

#### What it does mechanically — uncoupled block-diagonal assembly

The section order is `baseOrder + numMats` ([SectionAggregator.cpp:660-669](#)).
The response ordering is **base-section codes first, then the aggregated codes in
the order given** ([SectionAggregator.cpp:637-658](#), `getType()`).

**Stiffness** (`getSectionTangent`, [SectionAggregator.cpp:475-500](#)) copies the
base section's full (possibly coupled) tangent into the top-left block, then writes
**each addition's scalar tangent onto the diagonal only**:

```
(*ks)(i,i) = theAdditions[i-baseOrder]->getTangent();   // line 497
```

So the aggregated DOFs are **purely diagonal**: they do **not** couple with each
other, nor with the base section. Confirmed block-diagonal: the only off-diagonal
terms come from inside the base section's own tangent block; the additions never
write off-diagonal entries ([SectionAggregator.cpp:489-497](#)). `getStressResultant`
([SectionAggregator.cpp:593-614](#)) mirrors this — base block from
`theSection->getStressResultant()`, then each addition's scalar `getStress()`.
`getSectionFlexibility` inverts each diagonal addition tangent independently
(`1/k`, with a `1e14` singular guard, [SectionAggregator.cpp:550-559](#)).

> **Consequence:** SectionAggregator gives you a *block-diagonal* section. A fiber
> base section couples P-Mz-My among itself; the aggregated Vy/Vz/T springs are
> independent of everything. There is **no shear-flexure or shear-axial
> interaction** through the aggregator — if you need that, it must live inside a
> single coupled section/material, not across the aggregator seam.

#### Strain-rate / damping behavior (matters for the hinge use case)

`setTrialSectionDeformation()` ([SectionAggregator.cpp:427-450](#)) splits the
incoming deformation vector: the first `baseOrder` entries go to the base section,
the rest to the additions via

```
theAdditions[i-baseOrder]->setTrialStrain(def(i));   // line 447
```

This is the **single-argument** `setTrialStrain` — the `UniaxialMaterial` base
declares `setTrialStrain(double strain, double strainRate = 0)`
([UniaxialMaterial.h:59](#)), so the aggregated materials always receive
**strainRate = 0**. SectionAggregator therefore **does not pass strain rate** to
its uniaxial materials.

> **Net effect (combined with §2.2):** a `Viscous` / rate-dependent uniaxial
> material aggregated into a section is **inert as a dashpot on two counts**:
> (1) SectionAggregator never feeds it a rate (it calls the 1-arg
> `setTrialStrain`), and (2) even if it did, `ZeroLengthSection::update()` only
> calls `setTrialSectionDeformation(*v)` with deformations and **never passes a
> section strain rate** ([ZeroLengthSection.cpp:289](#)). So **do not** try to put
> a viscous spring in a `section Aggregator` for a `ZeroLengthSection` — it will
> produce zero velocity-proportional force. Velocity damping for a
> ZeroLengthSection hinge must come through Rayleigh (`-doRayleigh`, default ON
> for that element). Rate-independent uniaxial materials (Steel02, Hysteretic,
> Elastic, IMK*, Pinching, etc.) are unaffected and work exactly as intended.

#### State / commit, responses, gotchas

- **Commit cycle** ([SectionAggregator.cpp:671-717](#)): `commitState`,
  `revertToLastCommit`, `revertToStart` simply forward to the base section and
  every addition — all path-dependent state lives in the sub-materials/section,
  the aggregator holds none of its own (only scratch `e/s/ks/fs` work buffers).
- **Recorders** ([SectionAggregator.cpp:961-1016](#)): you can reach an aggregated
  material two ways — by code, `section ... -material <code> ...` /
  `addition <code> ...` (matches the first addition whose code equals the request,
  [SectionAggregator.cpp:985-990](#)), or by uniaxial **tag**,
  `material <matTag> ...` ([SectionAggregator.cpp:992-1000](#)). `section ...`
  forwards to the base section ([SectionAggregator.cpp:1004-1005](#)).
- **Order cap:** `maxOrder = 10` ([SectionAggregator.cpp:212](#)); base order +
  additions must not exceed it (hard `exit(-1)` otherwise,
  [SectionAggregator.cpp:264-268](#)).
- **Gotcha — duplicate / colliding codes are NOT checked.** Nothing prevents you
  from aggregating two materials on the same code (e.g. two `Mz`), or aggregating a
  code the base section already provides (e.g. base fiber section provides `P` and
  you also add a `P` spring). The result is **two independent entries with the same
  code in `getType()`**. Downstream consumers vary: `ZeroLengthSection`'s
  `setTransformation` maps *each* response row to its DOF, so two same-code rows
  both stamp the same element DOF — their stiffnesses **add** (often not what you
  want for a base+addition `P` collision), and code-based recorder lookups return
  only the **first** match. Order codes deliberately and avoid collisions with the
  base section.
- **Gotcha — ordering surprises.** Because the base section's codes come first and
  additions follow in input order, the section deformation/force vector layout is
  `[base codes... , addition codes in given order]`. If you write a custom
  element or recorder that assumes a fixed (P, Mz, My, Vy, Vz, T) order, it will be
  wrong — always consult `getType()`.

---

## 3. Damping & dashpots — consolidated

### 3.1 The two independent damping channels

For the material-based spring elements (ZeroLength, CoupledZeroLength,
TwoNodeLink) there are **two distinct ways to get damping**, and they must not be
confused:

1. **A velocity-dependent uniaxial material** (e.g. `Viscous`, `ViscousDamper`).
   This is a genuine dashpot: `F = C·v^α`. It works because `update()` feeds the
   material a *strain rate* via the 2-argument `setTrialStrain(strain, rate)`
   ([ZeroLength.cpp:819](#), [TwoNodeLink.cpp:663](#),
   [CoupledZeroLength.cpp:435](#)). The damping force comes out through the normal
   `getStress()` / `getResistingForce()` path — it is **always active**, no flag
   needed.

2. **Rayleigh `αM + βK` damping**, assembled in `getDamp()` from the element's
   stiffness/mass. This is gated by the `-doRayleigh` flag.

ZeroLengthSection and ZeroLengthND only ever offer channel (2) (Section) or
nothing (ND) — they have **no rate channel**, so a `Viscous` material is not even
an option there.

### 3.2 The `-doRayleigh` trap — per element

| Element | Flag | Default | What "off" means |
|---|---|---|---|
| ZeroLength | `-doRayleigh $f` | **0 (OFF)** ([ZeroLength.cpp:171](#)) | `getDamp()` uses material `getDampTangent()` (≈0 for elastic/plastic mats) → **stiffness-proportional Rayleigh silently does nothing** |
| ZeroLengthSection | `-doRayleigh $f` | **1 (ON)** ([ZeroLengthSection.cpp:76](#)) | If you set 0, `getDamp()` returns zero — *no* damping |
| ZeroLengthND | (none) | — | `getDamp()` always zero |
| CoupledZeroLength | `<useRayleigh>` | **0 (OFF)** ([CoupledZeroLength.cpp:81](#)) | same trap as ZeroLength |
| TwoNodeLink | `-doRayleigh` | **0 (OFF)** ([TwoNodeLink.cpp:149](#)) | material damping STILL applied; only the extra Rayleigh αM+βK is off |

> **Practical rule.** If you want global Rayleigh damping (`rayleigh αM βK ...`)
> to actually reach a ZeroLength / CoupledZeroLength spring, you **must** add
> `-doRayleigh 1`. If you forget, the spring contributes no βK damping and you
> get a silently under-damped model — a classic source of spurious high-frequency
> ringing across lumped-spring connections. ZeroLengthSection is the inverse:
> Rayleigh is on by default, which can *over*-damp a hinge you intended to be
> conservative.

### 3.3 Dashpot / absorbing-boundary recipe

To build a viscous dashpot (Lysmer–Kuhlemeyer style), use a **`Viscous`
uniaxialMaterial inside a ZeroLength or TwoNodeLink** on the desired direction.
The dashpot constant `C` (and exponent α) come from the material; the element
just projects it. No `-doRayleigh` needed. Example coefficients for a 1D
absorbing boundary: `C = ρ·c_p·A` (normal) and `ρ·c_s·A` (shear), with `c_p`,
`c_s` the P/S wave speeds.

---

## 4. Capability assessment for the four use cases

### 4.1 Contact (gap / normal / friction, unilateral)

**Use:** `ZeroLengthContact2D/3D` or `ZeroLengthContactASDimplex` for node-to-node
penalty contact; for soft/configurable gap behavior, a **ZeroLength with an
`ENT` (elastic-no-tension) or `Gap`/`ElasticPPGap` uniaxialMaterial** on the
normal direction plus a friction material on the tangent.

- **Stiffness:** the contact elements give a state-dependent tangent (open → 0,
  closed → Kn). A ZeroLength+ENT gives the same unilateral normal stiffness but
  with the friction coupling left to you (independent tangent springs, *no*
  normal-force-dependent friction cap).
- **Damping:** ZeroLengthContact* has no general dashpot; ASDimplex/Impact3D have
  internal damping. A ZeroLength gap spring can carry a `Viscous` material in
  parallel for contact damping.
- **Limitation / gap:** true Coulomb coupling (slip limit = μ·N) only exists in
  the dedicated contact elements, *not* in a hand-built ZeroLength (its springs
  are uncoupled). For **explicit** dynamics use `ZeroLengthContactASDimplex`
  (IMPL-EX makes the friction rate-independent and robust). This fork's
  Ladruno-contact plan extends exactly this element.

### 4.2 Phenomenological hinges (lumped plasticity)

**Use:** `ZeroLength` (independent per-DOF springs) for the common "rotational
spring + rigid translations" idealization; `ZeroLengthSection` when you need
**axial–moment coupling** from a fiber section; `ZeroLengthND` for fully coupled
multi-axial constitutive hinges; `TwoNodeLink` when the hinge must also transmit
P-Δ / shear-moment-arm effects across a finite offset.

- **Stiffness:** ZeroLength assembles `Σ t^T E t` per spring — perfect for a
  moment-rotation backbone (e.g. `Steel02`, `Hysteretic`, `IMKBilin`/`Pinching`
  on dir 6 = `Mz`). Put very large elastic springs on the DOFs that should be
  "rigid". ZeroLengthSection couples P and M through the section tangent.
- **Damping:** default `-doRayleigh 0` on ZeroLength means a stiffness-Rayleigh
  hinge contributes no damping — usually *fine* for a hinge (you don't want
  numerical damping smearing the hysteresis), but be deliberate. If you do want
  it, add `-doRayleigh 1`. ZeroLengthSection's default-on Rayleigh can
  unintentionally damp the hinge.
- **Limitation:** for the rigid DOFs you rely on a stiff penalty spring (see
  §4.4 caveats — conditioning). ZeroLength's springs are uncoupled, so a true
  P–M–V interaction surface needs ZeroLengthSection or ZeroLengthND.
- **Beam-column joints are a special case** — when the "hinge" is actually a
  panel zone (coupled panel-shear + bar-slip + interface, finite joint size),
  use a dedicated **joint element** (§7), not a single rotational ZeroLength.

#### Canonical full-6-DOF uncoupled hinge: `section Aggregator` + `ZeroLengthSection`

The standard pattern for a phenomenological hinge that needs **all 6 section DOFs
defined but only some of them coupled** is:

1. Build (or reuse) a **fiber base section** that couples `P-Mz-My` (the fiber
   integration gives the real axial-flexure interaction).
2. **Aggregate uncoupled `Vy`, `Vz`, `T` uniaxial springs** onto it with
   `section Aggregator ... -section $fiberTag` — the shear and torsion responses
   the fiber section does not provide ([SectionAggregator.h:32-35](#) — this is
   literally the design intent).
3. Feed the aggregated section to **one `ZeroLengthSection`**.

Result: a single zero-length hinge with a **6-DOF section** whose P-M block is
fiber-coupled and whose V/T springs are independent (block-diagonal, §2.7). This
is the canonical "fiber hinge with elastic shear/torsion" model.

```tcl
# 1) fiber base section couples P, Mz, My (tag 100)
section Fiber 100 -GJ 1.0e10 {
    # ... patch / layer / fiber definitions ...
}
# 2) elastic (rate-independent!) shear + torsion springs
uniaxialMaterial Elastic 201 [expr $G*$Avy]   ;# Vy
uniaxialMaterial Elastic 202 [expr $G*$Avz]   ;# Vz
uniaxialMaterial Elastic 203 [expr $G*$J]     ;# T
# 3) aggregate them onto the fiber section -> full 6-DOF uncoupled-shear section (tag 7)
section Aggregator 7 201 Vy 202 Vz 203 T -section 100
# 4) consume in a ZeroLengthSection hinge (i->j coincident nodes)
element zeroLengthSection 11 100 101 7      ;# -doRayleigh defaults ON
```

```python
ops.section('Fiber', 100, '-GJ', 1.0e10)   # + fiber/patch/layer calls
ops.uniaxialMaterial('Elastic', 201, G*Avy)
ops.uniaxialMaterial('Elastic', 202, G*Avz)
ops.uniaxialMaterial('Elastic', 203, G*J)
ops.section('Aggregator', 7, 201, 'Vy', 202, 'Vz', 203, 'T', '-section', 100)
ops.element('zeroLengthSection', 11, 100, 101, 7)
```

- **Section response order** is `[P, Mz, My  (from fiber) , Vy, Vz, T  (additions)]`
  — base codes first, additions in input order (§2.7,
  [SectionAggregator.cpp:637-658](#)). `ZeroLengthSection::setTransformation`
  maps each to its element DOF.
- **Do NOT use a `Viscous` material for the shear/torsion springs** here: neither
  SectionAggregator nor ZeroLengthSection passes strain rate (§2.7), so it would
  give zero damping force. Use `Elastic` (or a nonlinear rate-independent material
  if you want inelastic shear/torsion). Velocity damping for this hinge comes only
  from Rayleigh (`-doRayleigh`, default ON for ZeroLengthSection).
- **Pure rotational-spring hinges** (no axial-flexure coupling needed) are simpler
  with plain `ZeroLength` + a moment-rotation material on dir 6; reach for the
  Aggregator+ZeroLengthSection pattern when you specifically want the fiber P-M
  coupling plus a full 6-DOF section.

### 4.3 Absorbing boundaries (dashpots / Lysmer–Kuhlemeyer)

**Use:** `ZeroLength` (or `TwoNodeLink` if you want lumped mass too) with a
`Viscous` uniaxialMaterial per direction; tie the far node to a fixed support.

- **Damping:** this is the canonical dashpot path — `setTrialStrain(strain,rate)`
  ([ZeroLength.cpp:819](#)) gives the material the velocity it needs.
  **No `-doRayleigh` required** (the force is a stress, not a Rayleigh term).
- **Stiffness:** a pure dashpot has zero tangent stiffness (a `Viscous` material
  returns 0 stiffness, nonzero damp tangent) — fine for an explicit analysis;
  for implicit, the damping tangent enters via `getDamp()`.
- **Limitation / gap:** the ZeroLength+Viscous construction is the *manual*
  Lysmer dashpot — the lightweight, per-node approach most SSI/DRM models use.
  OpenSees does, however, ship **dedicated, engineered absorbing-boundary
  elements outside the zeroLength/link family** (not covered in depth here):
  `ASDAbsorbingBoundary2D/3D` (Lysmer–Kuhlemeyer with the free-field correction),
  `LysmerTriangle` + `LysmerVelocityLoader` (dashpot boundary element + input
  motion), and a full **PML** family (`PML2D/3D`, `PML2D_3/5/12`, `PML*VISCOUS`).
  Prefer those when you need a directionally-consistent, free-field-correct
  boundary; reach for ZeroLength+Viscous when you want a simple, tunable dashpot
  you control element-by-element. ZeroLengthSection / ZeroLengthND **cannot** do
  any of this (no rate channel). For the manual boundary, size each `Viscous` C by
  `ρ·c·A` and orient with `-orient`. (This is squarely in scope for the Ladruno
  explicit/SSI work — the dashpot path is verified above to feed strain rate
  correctly.)

### 4.4 Out-of-OpenSees typed constraints (penalty / MPC emulation, rigid links)

**Use:** `ZeroLength` with a very stiff `Elastic` uniaxialMaterial on the
constrained DOF(s) to emulate `equalDOF` / a rigid link as a **penalty
constraint**.

- **How:** `element zeroLength $t $i $j -mat $kBig -dir $d`. The relative
  deformation `u_j − u_i` along dir `d` is penalized by `kBig`. Multiple `-mat
  -dir` pairs lock multiple DOFs. For a rigid offset/moment arm use
  **TwoNodeLink** with `-shearDist` so the constraint reaction lands at the right
  point.
- **Stiffness considerations:** `kBig` must be large vs. the surrounding
  structural stiffness but not so large it wrecks conditioning (rule of thumb
  10^3–10^6 × neighbor stiffness). Unlike a true transformation/Lagrange
  constraint, the penalty spring leaves a small residual relative displacement
  and adds a stiff mode to the system.
- **Damping considerations:** a stiff penalty spring with global Rayleigh βK and
  **`-doRayleigh 0` contributes no βK damping** (the §3.2 trap) — which is
  actually desirable here, because you do *not* want a huge βK·kBig term
  injecting artificial damping/forces. Leave `-doRayleigh` at its default 0.
- **Limitation / gap:** penalty springs are an *approximation* of a constraint —
  for exact kinematic constraints prefer the real constraint handlers
  (`equalDOF`, `rigidLink`, `Transformation`/`Lagrange` constraint handler). The
  spring approach is the right tool when you need a *typed, tunable* (slightly
  flexible, or nonlinear) connection that constraint handlers can't express, or
  when you want the connection to also report a force.

---

## 5. Spring-defining materials — the constitutive half

A `ZeroLength` / `TwoNodeLink` is only a **frame**: it projects a relative
deformation (and, on the rate-capable elements, a relative velocity) onto a
`uniaxialMaterial` and sums `tᵀ E t`. **The material is what makes the spring
behave like contact, a dashpot, or a hinge.** All paths below are under
`SRC/material/uniaxial/`.

The decisive column is **rate-dependent?** — i.e. does the material actually
consume the `strainRate` arg of `setTrialStrain(strain, rate)` and return a
velocity-dependent force. Only rate-dependent materials act as true dashpots, and
only inside the rate-capable elements (ZeroLength, CoupledZeroLength, TwoNodeLink
— **not** ZeroLengthSection/ND, §3).

### 5.1 Contact / gap (unilateral) materials

| Material | Command | Law (1-line) | Rate-dep? |
|---|---|---|---|
| **ElasticPPGap** | `uniaxialMaterial ElasticPPGap $tag $E $Fy $gap <$eta $damage>` | Elastic-perfectly-plastic with an initial gap; zero force until \|ε\|>gap, then `E` to `Fy`, then `eta·E`. Sign of Fy/gap sets tension- vs compression-only. `getInitialTangent`=0 when gap open ([EPPGapMaterial.cpp:147-200](#)). | **No** |
| **ENT** | `uniaxialMaterial ENT $tag $E <$a $b>` | Elastic-no-tension: `σ=E·ε` in compression, 0 (or soft `a·E·tanh(b·ε)`) in tension. Switch at the origin, no gap distance ([ENTMaterial.cpp:110-117](#)). | **No** |
| **HookGap** | `uniaxialMaterial HookGap $tag $E $gap` or `$tag $E $gapN $gapP` | Double-sided elastic gap: zero force inside `[gapN,gapP]`, elastic `E` beyond either edge ([HookGap.cpp:139-164](#)). | **No** |
| **ImpactMaterial** | `uniaxialMaterial ImpactMaterial $tag $K1 $K2 $Delta_y $gap` | Compression-only bilinear **hysteretic** contact (K1 then K2 past Δy; separate unload path). Dissipative but the dissipation is hysteretic, **not viscous** ([ImpactMaterial.cpp:131-157](#)). | **No** |
| **HyperbolicGapMaterial** | `uniaxialMaterial HyperbolicGapMaterial $tag $Kmax $Kur $Rf $Fult $gap` | Compression-only hyperbolic (Duncan-Chang) backbone asymptoting to `Fult`, unload/reload at `Kur` ([HyperbolicGapMaterial.cpp:134-155](#)). | **No** |
| **Hertzdamp** | `uniaxialMaterial Hertzdamp $tag $Kh $xiNorm $gap <$n>` | Nonlinear Hertz contact `F=(Kh+ξv)·d^n` with velocity damping **only on approach** ([Hertzdamp.cpp:162-172](#)). | **YES** |
| **ViscoelasticGap** | `uniaxialMaterial ViscoelasticGap $tag $K $C $gap` | Kelvin-Voigt contact `F=K·d + C·v` while closed ([ViscoelasticGap.cpp:157](#)). | **YES** |
| **JankowskiImpact** | `uniaxialMaterial JankowskiImpact $tag $Kh $xi $mEff $gap <$n>` | Nonlinear viscoelastic (Hertz-damp/Jankowski) impact, damping only during approach ([JankowskiImpact.cpp:165-167](#)). | **YES** |

`Kh`/`Kmax`/`K`/`K1`/`E` is the **contact penalty stiffness** — set it large vs.
neighbouring members to approach impenetrability; `gap` is the standoff (negative
for the compression-only conventions of Impact/EPPGap). The three viscoelastic
gaps (Hertzdamp, ViscoelasticGap, JankowskiImpact) are the only contact materials
that add **viscous** damping; everything else is elastic or hysteretic.

### 5.2 Dashpot / viscous (rate-dependent) materials

| Material | Command | Model | `getTangent` / `getDampTangent` | Rate channel |
|---|---|---|---|---|
| **Viscous** | `uniaxialMaterial Viscous $tag $C $alpha <$minVel>` | Pure dashpot `F=C·sgn(v)·\|v\|^α` ([ViscousMaterial.cpp:122-138](#)) | tangent **0**, dampTangent **≠0** = `αC\|v\|^(α-1)` ([ViscousMaterial.cpp:141-162](#)) | uses `strainRate` |
| **ViscousDamper** | `uniaxialMaterial ViscousDamper $tag $K $C $Alpha <$LGap …>` | **Maxwell** (linear spring K series nonlinear dashpot), internal adaptive ODE solver ([ViscousDamper.cpp:173-205](#)) | tangent **0**, dampTangent **0**; force via `getStress()` | `strainRate` + `ops_Dt` |
| **Maxwell** | `uniaxialMaterial Maxwell $tag $K $C $Alpha Length $L` | Maxwell viscoelastic, closed-form exponential relaxation `tR=(C/Lᵅ)/K` ([Maxwell.cpp:142-155](#)) | tangent **K** (≠0), dampTangent 0 | `dStrain`/`ops_Dt` (NOT the `strainRate` arg) |

> **Two surprises worth internalizing:**
> 1. **`Viscous` has zero static stiffness.** `getTangent()`/`getInitialTangent()`
>    return 0; the damping enters only through `getDampTangent()`. Used **alone** in
>    a ZeroLength it makes the static tangent singular — **parallel it with an
>    elastic spring** (a second `-mat`/`-dir` pair, a `Parallel` material, or rely
>    on another element to supply stiffness on that DOF).
> 2. **`ViscousDamper` and `Maxwell` are rate-dependent but report
>    `getDampTangent()==0`** and integrate off the global `ops_Dt` rather than the
>    `strainRate` argument. They produce the right force **only in a transient
>    analysis with a defined time step**, and their rate effect is *not* exposed
>    through the damp-tangent channel the way `Viscous`'s is. For a clean,
>    tangent-visible dashpot (e.g. an absorbing boundary), prefer **`Viscous`**.

### 5.3 Phenomenological-hinge (degrading / pinching) materials — all rate-INDEPENDENT

| Material | Command gist | Captures |
|---|---|---|
| **Steel01 / Steel02** | `Steel01 $tag $Fy $E0 $b …` / `Steel02 $tag $Fy $E0 $b $R0 $cR1 $cR2 …` | Bilinear kinematic (01) / Giuffré-Menegotto-Pinto smooth bilinear with Bauschinger (02). |
| **Hysteretic** | `Hysteretic $tag $m1p $r1p $m2p $r2p <$m3p $r3p> $m1n … $pinchX $pinchY $damage1 $damage2 <$beta>` | 3-point multilinear backbone + pinching + ductility/energy strength & stiffness degradation ([HystereticMaterial.cpp:54-58](#)). |
| **Pinching4** | `Pinching4 $tag <4-pt pos> <4-pt neg> $rDispP $fForceP $uForceP … $gK1.. $gD1.. $gF1.. $gE $dmgType` | 4-point pinched backbone with **three independent** cyclic degradation laws (unload-stiffness, reload-stiffness, strength) ([Pinching4Material.cpp:88-113](#)). Richest model here. |
| **Bilin / ModIMKPeakOriented(02) / ModIMKPinching(02)** | `Bilin $tag $Ke $AsPos $AsNeg $My_pos $My_neg $LamS $LamC $LamA $LamK $cS $cC $cA $cK $thetaP… ` | Ibarra-Medina-Krawinkler hinge (bilinear / peak-oriented / pinched) with the **four IMK deterioration modes** (basic-strength S, post-cap C, accelerated A, unloading-stiffness K) ([Bilin.cpp:75](#), [ModIMKPeakOriented.cpp:75](#), [ModIMKPinching.cpp:73](#)). |
| **IMKBilin / IMKPeakOriented / IMKPinching** | `IMK<X> $tag $Ke …` | Newer reformulation of the IMK family (same three flavours) ([IMKBilin.cpp:62](#), [IMKPeakOriented.cpp:63](#), [IMKPinching.cpp:63](#)). |
| **SelfCentering** | `SelfCentering $tag $k1 $k2 $ActF $beta <$SlipDef $BearDef $rBear>` | Flag-shaped recentering hysteresis (SMA / post-tensioned) ([SelfCenteringMaterial.cpp:48-70](#)). |
| **HystereticSM** | `HystereticSM $tag …` | Extended multilinear hysteretic (more backbone points + degradation options than `Hysteretic`). |

> **All of §5.3 is rate-independent** — none reads `strainRate`. Damping for these
> hinges must come from **Rayleigh** (element/nodal/modal), never the material.
> **Quirk:** the `ModIMK*` materials override `getDampTangent()` to return the
> *current stiffness tangent* `ekP`, not a velocity derivative
> ([ModIMKPeakOriented.cpp:744-747](#)). If an element queries the damp tangent it
> gets the elastic-branch stiffness, not zero — harmless for ZeroLength's default
> path but worth knowing.

---

## 6. Bearing / isolator elements — purpose-built coupled springs

Under `SRC/element/elastomericBearing/` and `SRC/element/frictionBearing/`. These
are the **richer cousins of TwoNodeLink**: two-node, 6-DOF/node spring elements
whose **defining feature is a bidirectionally coupled shear law**.

### 6.1 The defining pattern: coupled shear, uncoupled everything else

- **The two shear directions are COUPLED** through one yield surface / hysteretic
  law (plasticity, Bouc-Wen, or friction-pendulum). This circular interaction —
  yielding/sliding in one horizontal direction reduces capacity in the orthogonal
  one — is the whole reason the family exists. A generic spring assembly cannot
  enforce it. (Plasticity = 2D radial return on the shear-force vector,
  [ElastomericBearingPlasticity3d.cpp:479-511](#); Bouc-Wen = 2-vector `z` solved
  by a 2×2 Newton, [ElastomericBearingBoucWen3d.cpp:505-577](#).)
- **Axial (P), torsion (T), and the two moments (My, Mz) are UNCOUPLED**, each
  handled by a separately-assigned `uniaxialMaterial` via `-P/-T/-My/-Mz`
  ([ElastomericBearingBoucWen3d.cpp:285-289](#)). The 3D commands require **4**
  materials, the 2D commands **2** ([TclElastomericBearingPlasticityCommand.cpp:163,395](#)).
  These materials **do** receive strain rate (`setTrialStrain(ub, ubdot)`,
  [ElastomericBearingBoucWen3d.cpp:493](#)), so rate-dependent uniaxials work on
  the uncoupled DOFs; the coupled shear law does its own rate handling.

### 6.2 Elastomeric (rubber) bearings

| Element | Command gist | Coupled shear law | `-doRayleigh` |
|---|---|---|---|
| **ElastomericBearingPlasticity 2d/3d** | `elastomericBearing $tag $i $j $kInit $qd $alpha1 $alpha2 $mu -P $m <-T -My> -Mz $m <-orient…> <-shearDist $r> <-doRayleigh> <-mass $m>` | Bilinear, 2D radial-return on shear force | **0 (OFF)** |
| **ElastomericBearingBoucWen 2d/3d** (+ **Mod3d**) | adds `$eta $beta $gamma` + `<-iter $maxIter $tol>` (default 25 / 1e-12) | Smooth Bouc-Wen, 2-vector z | **0 (OFF)** |
| **ElastomericBearingUFRP2d** | `… $uy $a1..$a5 $b $c $eta $beta $gamma -P $m -Mz $m …` | Bouc-Wen-type, 7-term polynomial backbone (U-shaped FRP) | 0 |
| **ElastomericX / LeadRubberX / HDR** (X-series) | **property-based**: rubber `Gr,Kbulk`, inner/outer dia, shim & rubber thickness, layer count → element derives stiffness/buckling/cavitation internally ([ElastomericX.cpp:265-304](#)) | internal coupled Bouc-Wen (+ lead-core heating for LeadRubberX) | **no flag** — explicit viscous term `cd` instead ([ElastomericX.cpp:660-663](#)) |

The **X-series** is the physical-property branch (with `.tcl` helpers): instead of
`kInit/qd/alpha` you give geometry and material moduli, and it models
vertical↔horizontal coupling (axial load modulates horizontal stiffness) that the
basic elements do not.

### 6.3 Friction (sliding) bearings

| Element | Command gist | Coupled shear law | `-doRayleigh` |
|---|---|---|---|
| **FlatSliderSimple 2d/3d** | `flatSliderBearing $tag $i $j $frnMdlTag $kInit -P $m <-T -My> -Mz $m <-orient…> <-shearDist $r> <-doRayleigh> <-iter …>` | Coulomb sliding via a `FrictionModel`, circular slip surface, **no restoring** | **0 (OFF)**; `shearDist` default **0.0** (slider at bottom node) |
| **SingleFPSimple 2d/3d** | adds `$Reff` before `$kInit`; `<-inclVertDisp>` | Single friction pendulum: friction + pendulum restoring of radius `Reff` (`R=√(Reff²−u²)`, [SingleFPSimple3d.cpp:523-536](#)) | **0 (OFF)** |
| **TripleFrictionPendulum / …X / TFP_Bearing(2d) / MultiFP2d / TPB1D** | `… $frn1 $frn2 $frn3 -P -Mz -My -T … $L1 $L2 $L3 $Ubar1..3 $W $uy …` | Three series pendulum surfaces, each its own `FrictionModel`, multi-stage backbone | Newton-iterated |
| **RJWatsonEQS 2d/3d** | `… $frnMdlTag $kInit -P $m -Vy $m <-Vz -T -My> -Mz $m …` | Flat slider + **mass-energy-regulator restoring springs** on `-Vy/-Vz` | **0 (OFF)** |
| **FPBearingPTV** | property/physical FP | Single-FP with μ varying by **Pressure, Temperature (frictional heating), Velocity** | — |

**Friction models** (`frictionModel/`, fed live normal force + slip velocity each
step, [FrictionModel.h:54-58](#)): `Coulomb` (constant μ), `VelDependent`
(`μ=μF−(μF−μS)e^(−rate·|v|)`), `VelDepMultiLinear` (piecewise μ(v)), `VelPressureDep`
(Constantinou μ(v,p)), `VelNormalFrcDep` (μ(v,N)).

### 6.4 When to use a bearing element vs. generic CoupledZeroLength + materials

Reach for a purpose-built bearing when you need any of:

1. **A true bidirectional shear yield surface** (circular interaction). `CoupledZeroLength`
   couples only two chosen DOFs through *one* material and cannot host the
   friction/pendulum laws.
2. **Friction-pendulum geometry** — restoring stiffness that scales with axial load
   and instantaneous radius `√(Reff²−u²)`, optional vertical-drop coupling. No
   generic-spring equivalent.
3. **Velocity/pressure/temperature-dependent friction** via the pluggable
   `FrictionModel` — impossible from rate-independent uniaxial springs.
4. **Physical-property parameterization** with load-coupled stiffness, buckling,
   cavitation, lead-core heating (the X-series).
5. **Built-in bearing P-Δ / large-displacement kinematics** with the standard
   `-orient`/`-shearDist` slip-plane placement.

If you only need two *independent* nonlinear springs, or one coupled pair through
an ordinary hysteretic material, `twoNodeLink` / `zeroLength` / `CoupledZeroLength`
are lighter and sufficient. The bearing family pays off precisely when the shear
coupling law, the friction model, or the geometry-derived stiffness is the point.
Cross-cutting: `-doRayleigh` defaults **OFF** for the whole family; the X-series
has no `-doRayleigh` (it uses an explicit `cd` viscous term).

---

## 7. Joint elements — panel-zone spring bundles

Under `SRC/element/joint/`. A joint element is **not a single spring but a small
structural subassembly**: it models a beam-column panel zone as a *coupled system*
of panel-shear + bar-slip + interface springs, with the correct rigid-link
kinematics of a finite-size joint. This is the specialized cousin of the §4.2
phenomenological hinge — reach for it when the hinge *is* a joint.

### 7.1 Two modeling philosophies

**(A) Constraint-based extra-DOF joints — `Joint2D`, `Joint3D`** (Altoontash &
Deierlein 2004). The element **auto-creates an internal central node** with an
extra rotational DOF representing panel-zone shear distortion, then installs
`MP_Joint2D`/`MP_Joint3D` multi-point constraints between that central node and
each external node to enforce the rigid-link panel kinematics. The panel
hysteresis lives in a single `uniaxialMaterial` rotational spring on the extra DOF;
optional interface springs sit at each external connection. The element is
essentially a *spring bundle + constraint generator* — panel stiffness comes
entirely from the springs ([Joint2D.cpp:1103-1117](#)).

**(B) Monolithic multi-spring super-elements — `BeamColumnJoint2d/3d`,
`LehighJoint2d`, `Inno3DPnPJoint`** (Lowes & Altoontash 2003 / Mitra & Lowes
2007). A single 4-node finite-area element internally assembles **13 springs**
([BeamColumnJoint2d.h:30-34](#)): 8 bar-slip (anchorage) + 4 interface-shear + 1
shear-panel, all condensed internally via a compatibility matrix + static
condensation so the element exposes only the external nodal DOFs. **No internal
node, no MP constraints** — all coupling is inside the element
([BeamColumnJoint2d.cpp:1388-1396](#)).

The two outliers: `ElasticTubularJoint` (elastic offshore K/T/Y tubular joint,
geometric inputs, **no materials**) and `Inno3DPnPJoint` (a specific 32-spring 3D
plug-and-play steel connection).

### 7.2 Comparison

| Element | Nodes | Total DOF | Internal node? | # materials & what they model | Constraint-based? | Reference |
|---|---|---|---|---|---|---|
| **Joint2D** | 4 ext + 1 auto | 16 (4×3 + central 4) | **Yes** ([Joint2D.cpp:678](#)) | 1 panel-shear (mandatory) + up to 4 optional interface springs; **NULL interface mat ⇒ that end is fixed** ([Joint2D.cpp:688-700](#)) | **Yes** — 4× MP_Joint2D ([Joint2D.cpp:713-735](#)) | Altoontash 2004 |
| **Joint3D** | 6 ext + 1 auto | 45 (6×6 + central 9) | **Yes** ([Joint3D.cpp:321](#)) | 3 panel-shear springs (MatX/Y/Z, shear in 3 planes) ([Joint3D.cpp:332-348](#)) | **Yes** — 6× MP_Joint3D ([Joint3D.cpp:362-369](#)) | Altoontash 2004 (3D) |
| **BeamColumnJoint2d** | 4 | 12 (4×3) | No (condensed) | **13**: 8 bar-slip + 4 interface + 1 panel ([BeamColumnJoint2d.h:30-34](#)) | No (compat + condensation) | Lowes & Altoontash 2003 |
| **BeamColumnJoint3d** | 4 | 24 (4×6) | No | **13** (same triad) ([BeamColumnJoint3d.h:31](#)) | No | Lowes & Altoontash 2003 (3D) |
| **LehighJoint2d** | 4 | 12 (9 basic) | No | **9** springs ([LehighJoint2d.cpp:78-93](#)) | No (avp/apq compat) | Seo, Lehigh 2008 (reduced LA) |
| **ElasticTubularJoint** | 2 | 6 | No | **0** — geometric only (brace/chord D, t, angles, E) ([ElasticTubularJoint.cpp:71-128](#)) | No | Offshore tubular LJF (MSL/Fessler) |
| **Inno3DPnPJoint** | 5 | 30 | No | **32** springs ([Inno3DPnPJoint.cpp:784](#)) | No | 3D plug-and-play steel connection |

Joint2D/3D additionally accept optional per-spring `DamageModel`s (an
Altoontash-specific feature absent from the monolithic family,
[Joint2D.cpp:561-567](#)).

### 7.3 MP_Joint2D / MP_Joint3D (the constraint behind philosophy A)

These are **`MP_Constraint` subclasses, not elements**
([MP_Joint2D.h:48](#)). They impose the **large-deformation rigid-link
panel-zone constraint** between the auto-created central node (the *retained* node,
which carries the extra DOF) and each external node (the *constrained* node).
`MP_Joint2D` requires the retained node to have **4 DOF** and the constrained node
**3 DOF** ([MP_Joint2D.cpp:91](#)); the constraint is the rigid-arm relation
`u = u_c − Δy·θ`, `v = v_c + Δx·θ` ([MP_Joint2D.cpp:149-152](#)). Crucially it is
**time-varying and re-forms from current geometry** when `LargeDisplacement` is
set ([MP_Joint2D.cpp:127](#)) — giving the *correct geometric nonlinearity* of the
panel, which a plain `equalDOF`/`rigidLink` cannot. A `FixedEnd` flag chooses
whether the external rotational DOF is released or rigidly tied to the panel.

### 7.4 Damping — none of these are dashpots

> **There is no `-doRayleigh` flag anywhere in `SRC/element/joint/`** (grep returns
> nothing), and **none of these elements pass strain rate to their springs** —
> every `setTrialStrain` call is single-argument ([Joint2D.cpp:1015](#),
> [BeamColumnJoint2d.cpp:798](#), [LehighJoint2d.cpp:379](#),
> [Inno3DPnPJoint.cpp:1787](#)). These are **quasi-static / cyclic** models;
> viscous damping is not modeled inside the element. Any damping comes only from
> model-level Rayleigh via the base-class `getDamp()`. Use rate-independent
> hysteretic materials (Steel02, Hysteretic, Pinching4) for the springs — a
> `Viscous` material would be inert here, exactly as in ZeroLengthSection (§5.2).

### 7.5 When to use a joint element vs. a zeroLength rotational hinge

Use a plain **`zeroLength`/`zeroLengthSection` rotational spring** when you need a
single lumped moment-rotation hysteresis between two coincident nodes (a member
plastic hinge): one uncoupled DOF, no geometry.

Use a **joint element** when the panel zone is a *coupled system* — panel-zone
shear **plus** bar-slip/anchorage **plus** interface slip, over a *finite joint
area* with rigid-link kinematics as the frame deforms. A single rotational
`zeroLength` cannot represent (a) the finite physical size (rigid offsets that move
the beam/column ends with the panel), (b) the separate bar-slip and interface
mechanisms the Lowes-Altoontash 13-spring model isolates, or (c) the
large-deformation coupling between panel rotation and the corner nodes — enforced
exactly by the internal `MP_Joint` constraints (philosophy A) or the internal
compatibility + condensation (philosophy B).

### 7.6 Gotchas

- **Joint2D/3D auto-create the internal central node** — you supply a *new, unused*
  tag (`NodC`); the parser rejects it if that node already exists
  ([Joint2D.cpp:87-93](#)). Do **not** define it yourself.
- **External nodes must have exactly 3 DOF (Joint2D) / 6 DOF (Joint3D)** — checked
  in `setDomain` ([Joint2D.cpp:641-645](#)); the central node is created with 4 / 9
  DOF automatically.
- **Node order must form a parallelogram** in Joint2D (I-J-K-L around the panel),
  else "can not construct a parallelogram" error ([Joint2D.cpp:663-674](#)).
- **NULL interface material fixes that end** (`fixedEnd[i]=1`,
  [Joint2D.cpp:689-691](#)) — pass matTag `0` to rigidly connect an end. The
  **panel spring `MatC` is mandatory** (element `exit(-1)`s without it,
  [Joint2D.cpp:699](#)).
- **Two parser paths** coexist (legacy `TclModelBuilder_addJoint2D` and the newer
  `OPS_Joint2D`, [Joint2D.cpp:53](#)) with the same grammar.

#### Command signatures
```
element Joint2D $tag $NodI $NodJ $NodK $NodL $NodC $MatC $LrgDsp
element Joint2D $tag $NodI $NodJ $NodK $NodL $NodC $MatI $MatJ $MatK $MatL $MatC $LrgDsp <-damage ...>
element Joint3D $tag $NodI..$NodN $NodC $MatX $MatY $MatZ $LrgDsp <-damage $DmgX $DmgY $DmgZ>
element beamColumnJoint $tag $n1 $n2 $n3 $n4 $matTag1 ... $matTag13 <$Hgtfac $Wdtfac>
element LehighJoint $tag $n1 $n2 $n3 $n4 $matTag1 ... $matTag9
element ElasticTubularJoint $tag $iNode $jNode $braceD $braceAngle $E $chordD $chordT $chordAngle
element Inno3DPnPJoint $tag $Node1..$Node5 $Spring01..$Spring32
```

---

## 8. Quick reference / cheat-sheet

### Phenomenological rotational hinge (2D, dir 6 = Mz)
```tcl
# Tcl
element zeroLength 10 100 101 -mat 5 -dir 6 -doRayleigh 1
```
```python
# openseespy
ops.element('zeroLength', 10, 100, 101, '-mat', 5, '-dir', 6, '-doRayleigh', 1)
```

### Fiber lumped-plastic hinge with P-M coupling
```tcl
element zeroLengthSection 11 100 101 7 ;# section 7; -doRayleigh defaults ON
```

### section Aggregator — full 6-DOF uncoupled-shear hinge section
```tcl
# Elastic (rate-independent) Vy/Vz/T springs aggregated onto fiber section 100
uniaxialMaterial Elastic 201 [expr $G*$Avy]
uniaxialMaterial Elastic 202 [expr $G*$Avz]
uniaxialMaterial Elastic 203 [expr $G*$J]
section Aggregator 7 201 Vy 202 Vz 203 T -section 100
element zeroLengthSection 11 100 101 7
```
```python
ops.uniaxialMaterial('Elastic', 201, G*Avy)
ops.uniaxialMaterial('Elastic', 202, G*Avz)
ops.uniaxialMaterial('Elastic', 203, G*J)
ops.section('Aggregator', 7, 201, 'Vy', 202, 'Vz', 203, 'T', '-section', 100)
ops.element('zeroLengthSection', 11, 100, 101, 7)
```

### section Aggregator — standalone uncoupled springs (no base section)
```tcl
# P + Mz hinge section from two independent uniaxial materials (no fiber base)
section Aggregator 8 5 P 6 Mz       ;# mat 5 -> axial, mat 6 -> bending z
element zeroLengthSection 12 100 101 8
```
```python
ops.section('Aggregator', 8, 5, 'P', 6, 'Mz')
ops.element('zeroLengthSection', 12, 100, 101, 8)
```

### Beam-column joint panel zone (Joint2D, extra-DOF + auto central node)
```tcl
# 4 external nodes 1-2-3-4 form a parallelogram; node 99 is NEW (auto-created).
# mat 0 on an end = rigid; matC (10) = mandatory panel-shear hysteresis.
uniaxialMaterial Hysteretic 10 ... ;# panel-zone shear spring (rate-independent)
element Joint2D 70 1 2 3 4 99 0 0 0 0 10 1   ;# interface ends rigid, LrgDsp=1
```
```python
ops.element('Joint2D', 70, 1,2,3,4, 99, 0,0,0,0, 10, 1)
```

### Monolithic 13-spring joint (BeamColumnJoint, condensed, no internal node)
```tcl
element beamColumnJoint 71 1 2 3 4 \
    $m1 $m2 $m3 $m4 $m5 $m6 $m7 $m8 $m9 $m10 $m11 $m12 $m13
# m1-8 bar-slip, m9-12 interface-shear, m13 shear-panel (Lowes-Altoontash 2003)
```

### Viscous dashpot / absorbing boundary (3D normal, fixed far node)
```tcl
uniaxialMaterial Viscous 90 [expr $rho*$cp*$A] 1.0   ;# C, alpha
element zeroLength 20 200 201 -mat 90 -dir 1 -orient 1 0 0 0 1 0
# no -doRayleigh needed; force is velocity-proportional
```
```python
ops.uniaxialMaterial('Viscous', 90, rho*cp*A, 1.0)
ops.element('zeroLength', 20, 200, 201, '-mat', 90, '-dir', 1,
            '-orient', 1,0,0, 0,1,0)
```

### Node-to-node penalty contact + friction (2D)
```tcl
element zeroLengthContact2D 30 300 301 $Kn $Kt $fs -normal 1.0 0.0
```

### Explicit-friendly contact (IMPL-EX)
```tcl
element zeroLengthContactASDimplex 31 300 301 $Kn $Kt $mu -orient 1 0 0
```

### Penalty "equalDOF" (lock dir 1 & 2 between coincident nodes)
```tcl
uniaxialMaterial Elastic 99 1.0e9
element zeroLength 40 400 401 -mat 99 99 -dir 1 2
```

### Material-built contact spring (gap normal + uncoupled friction tangent)
```tcl
uniaxialMaterial ElasticPPGap 70 $Kn 0.0 -0.001   ;# compression-only gap, standoff 1mm
uniaxialMaterial Steel01     71 $fs $Kt 0.0       ;# crude friction (NOT coupled to N)
element zeroLength 32 320 321 -mat 70 71 -dir 1 2 -orient 1 0 0 0 1 0
# NB: slip limit is fixed fs, not mu*N. For true Coulomb coupling use ZeroLengthContact*.
```

### Viscous dashpot paralleled with an elastic spring (avoids singular tangent)
```tcl
uniaxialMaterial Viscous 80 $C 1.0     ;# pure dashpot: ZERO static stiffness
uniaxialMaterial Elastic 81 $kStab     ;# small stabilizing stiffness in parallel
element zeroLength 33 330 331 -mat 80 81 -dir 1 1   ;# both on dir 1 -> stiffnesses add
```

### Friction-pendulum isolator (single, 3D)
```tcl
frictionModel VelDependent 1 $muSlow $muFast $rate
uniaxialMaterial Elastic 90 $kV        ;# vertical (P)
uniaxialMaterial Elastic 91 1.0e10     ;# torsion / moments (rigid)
element singleFPBearing 60 600 601 1 $Reff $kInit \
    -P 90 -T 91 -My 91 -Mz 91 -orient 0 0 1 1 0 0 -shearDist 0.0
```

### Elastomeric (lead-rubber) bearing, bilinear coupled shear (3D)
```tcl
uniaxialMaterial Elastic 95 $kV
uniaxialMaterial Elastic 96 1.0e10
element elastomericBearing 61 610 611 $kInit $qd $alpha1 0.0 0.0 \
    -P 95 -T 96 -My 96 -Mz 96 -orient 0 0 1 1 0 0 -shearDist 0.5
```

### Base isolator / damper with mass and P-Delta (TwoNodeLink)
```tcl
element twoNodeLink 50 500 501 -mat 60 61 62 -dir 1 2 3 \
    -pDelta 0.5 0.5 0.5 0.5 -shearDist 0.5 0.5 -mass 2.0 -doRayleigh
```
```python
ops.element('twoNodeLink', 50, 500, 501, '-mat', 60,61,62, '-dir', 1,2,3,
            '-pDelta', 0.5,0.5,0.5,0.5, '-shearDist', 0.5,0.5,
            '-mass', 2.0, '-doRayleigh')
```

---

## Appendix — direction code (all material-based elements)

`-dir`/`-dof` value (1-based input, stored 0-based):

| input | local DOF | meaning |
|---|---|---|
| 1 | dx | axial / translation along local x |
| 2 | dy | translation along local y (in 2D-3dof, "2" is remapped to Mz in ZeroLength) |
| 3 | dz | translation along local z |
| 4 | rx | rotation about local x (torsion) |
| 5 | ry | rotation about local y |
| 6 | rz | rotation about local z |

Local axes are set by `-orient` (default global x, y); with coincident nodes the
orientation is *only* what you supply. TwoNodeLink additionally derives local-x
from the node geometry when the nodes are not coincident.
