---
title: ADR — Concentrated-plasticity IMK beam macro element (LadrunoIMKBeam)
project: Ladruno
status: proposed
priority: medium
owner: nmora
tags:
  - implementation
  - element
  - beam-column
  - concentrated-plasticity
  - imk
  - adr
---

# ADR — Concentrated-plasticity IMK beam (`LadrunoIMKBeam`)

**Status:** proposed (2026-06-02) · **Registry:** `Element`, classTag **33003**
(next free Element tag after BezierTri6=33000, BezierTet10=33001,
LadrunoBrick=33002; class-tag bands are per-registry) · **Consumes:** any two
3D nodes (ndf=6), an existing `CrdTransf3d`, and pre-built `uniaxialMaterial`
moment-rotation laws.

> One-line verdict from the scoping discussion: this is a **usability +
> robustness** product, not novel mechanics. The uncoupled concentrated-plasticity
> model is already reachable in stock OpenSees via
> `forceBeamColumn` + `HingeRadau` + `section Aggregator{IMK on Mz, IMK on My,
> Elastic on P, Elastic on T}`. The macro earns its place by packaging that into
> one clean tag with: a **column-face offset**, **guaranteed-correct damping**,
> **single-tag recorders**, and **element-level softening robustness** — with
> **no n-factor, no artificial stiffness, no constraint-chaining** against the
> rigid diaphragm.

---

## 1. Context — why fiber/distributed hinges break under a rigid diaphragm

The motivating failure: modeling beams in a moment frame with a **rigid
diaphragm** plus `forceBeamColumn` + `HingeRadau` + a **fiber** hinge section.

- The diaphragm MPC enforces equal lateral displacement across the floor. Any
  axial-flexibility mismatch between framing members shows up as a **spurious
  axial force `P`** injected into the beams.
- A fiber hinge computes moment capacity *from* the axial-flexure interaction at
  the section. It dutifully shifts `Mz`/`My` capacity for that fake `P`, and when
  `P` approaches the section axial capacity the state determination misbehaves or
  fails to converge. The physics is contaminated **and** the solver is fragile.

The fix is the industry-standard ASCE-41 / Lignos-Krawinkler idiom:
**concentrated, uncoupled moment-rotation rotational springs.** Moment response
is a function of rotation only; the diaphragm's axial force flows through a
separate elastic `EA` path and never touches the moment law.

### 1.1 Why a dedicated element rather than the Python two-element recipe

The classic script (elastic beam + two `zeroLength` rotational springs) needs a
stack of workarounds because it glues **two independent OpenSees objects** in
series:

- the **n-factor** (Zareian-Medina 2010): inflate the spring to `n·6EI/L`,
  stiffen the elastic element to `(n+1)/n·I`, and modify the spring hardening to
  `α_s = α_mem/(1+n(1−α_mem))` so the *series* assembly recovers the true member
  stiffness;
- the resulting **artificial 10× stiffness** corrupts stiffness-proportional
  Rayleigh damping (the entire subject of Zareian-Medina 2010), forcing
  `-doRayleigh 0` discipline on the springs;
- the internal-node translational ties (`equalDOF`) **chain** against the
  `rigidDiaphragm` MPC, which the `Transformation` handler treats inconsistently.

A purpose-built element that owns the assembly **eliminates all three**: it adds
the hinge flexibility *directly* (exact series, no n-factor), so there is no
inflated stiffness (damping behaves) and no MPC (no chaining). See §2.

### 1.2 Explicit non-goals (v1)

- **No biaxial `My`-`Mz` coupling.** The two bending hinges are independent
  uniaxial laws (block-diagonal), by the same construction that buys axial
  immunity. Correct and standard for **beams** (strong-axis-dominated); **over**predicts
  capacity for **columns** under genuine biaxial demand. A coupled
  moment-interaction-surface hinge is a *different* element (`zeroLengthND` +
  coupled `nDMaterial`, or a future `LadrunoIMKColumn`) — recorded here as the
  follow-up, not a flag on this element.
- **No `P`-`M` interaction** (that is the whole point).
- **No material baked in.** Steel vs concrete is purely the choice of uniaxial
  tag (see §4).

---

## 2. Formulation — displacement-driven stiffness macro

Work in the element **basic system** (rigid-body modes removed; geometry +
offsets + P-Δ delegated to the `CrdTransf3d`). The 3D basic forces/deformations:

```
q = [ N , Mz_i , Mz_j , My_i , My_j , T ]ᵀ
v = [ Δ , θz_i , θz_j , θy_i , θy_j , φ ]ᵀ
```

The hinge is in **series** with the elastic interior — same force, deformations
add — so the basic flexibility adds directly:

```
F(q) = F_elastic + F_hinge(q)
```

with `F_hinge` contributing `1/k_IMK` **only** on the hinged bending diagonals.
No `n`-factor: the series is represented exactly.

### 2.1 Elastic interior (clear span Lc, Euler-Bernoulli v1)

```
axial    : f_NN = Lc/(E·A)
torsion  : f_TT = Lc/(G·J)
bending z: F_el,z = (Lc/(6·E·Iz)) · [[ 2, -1],[-1, 2]]      (θz_i, θz_j)
bending y: F_el,y = (Lc/(6·E·Iy)) · [[ 2, -1],[-1, 2]]      (θy_i, θy_j)
```

`Lc = L_full − offset_i − offset_j` (clear length, §3). Timoshenko shear
(`GA_s`) is a v1.x flag, adding `±Lc/(GA_s)`-scale terms to the bending blocks;
**default Euler-Bernoulli** (the standard IMK assumption).

### 2.2 Hinge flexibility and the per-axis internal solve

IMK is **displacement-driven** (`setTrialStrain(θ_h) → M, k_t`), so we solve for
the **hinge rotations** `θ_h` directly — no material inversion. Each bending axis
is an independent 2-equation system. For the `z` axis, unknowns
`θ_h = [θ_h,i, θ_h,j]`:

```
total compatibility :  v_z = F_el,z · q_z + θ_h
material law        :  q_z = [ M_i(θ_h,i) , M_j(θ_h,j) ]   (IMK)
```

Substituting gives the residual solved by a tiny Newton:

```
r_i = θ_h,i − v_z,i + F_el,z[i,i]·M_i(θ_h,i) + F_el,z[i,j]·M_j(θ_h,j)
r_j = θ_h,j − v_z,j + F_el,z[j,i]·M_i(θ_h,i) + F_el,z[j,j]·M_j(θ_h,j)

J = I + F_el,z · diag(k_t,i, k_t,j)         # 2×2, k_t = IMK tangent
Δθ_h = −J⁻¹ r                                # ~3 iters
```

Converged → basic-stiffness block for the axis:

```
K_b,z = (F_el,z + diag(1/k_t,i, 1/k_t,j))⁻¹
```

Axial and torsion are linear elastic (no hinge): `q_N = (EA/Lc)·Δ`,
`q_T = (GJ/Lc)·φ`. Assemble the 6×6 `K_b`, transform to the global 12×12 via the
`CrdTransf3d`.

**Robustness:** IMK softening (`k_t < 0`) is resolved *inside* the element's 2×2
solve, not dumped into the global tangent — a real advantage over the
two-element recipe where the negative stiffness reaches the assembler.

### 2.3 Hinge-present matrix

| `-hinge` | i end | j end |
|---|---|---|
| `both` | `F_hinge` added | `F_hinge` added |
| `i` | added | `θ_h,j ≡ 0` (rigid; M from elastic compatibility) |
| `j` | `θ_h,i ≡ 0` | added |

Per axis, `-matY` omitted ⇒ that axis is fully elastic (no hinge term, `My` from
`F_el,y` alone). Both axes absent ⇒ pure elastic beam `K_b = F_el⁻¹` (a useful
degenerate-mode self-check).

### 2.4 State, tangents, damping

- **Committed state:** the four scalar hinge rotations
  `θ_h ∈ {z_i, z_j, y_i, y_j}` plus each IMK material's own history
  (materials `commitState()` themselves). Trial reset on `revertToLastCommit`.
- **`getTangentStiff`:** §2.2 with current IMK tangents.
- **`getInitialStiff`:** §2.2 with each IMK initial tangent — **exact**, so
  initial-stiffness Rayleigh is clean. No inflated stiffness anywhere, so the
  Zareian-Medina βK pathology does not arise; document "prefer initial-stiffness
  Rayleigh" but it is no longer load-bearing.

---

## 3. Column-face offset — reuse `CrdTransf -jntOffset`

The hinge belongs at the **column face / RBS location**, not the joint
centerline. Rather than synthesize rigid stubs, reuse the existing rigid joint
offset in the geometric transformation:

- the element forms `F_elastic` on the **clear length** `Lc`;
- the `CrdTransf3d` rigid offset maps the clear-span ends (where the hinges sit)
  back to the real joint nodes `i`, `j`.

So the hinges land exactly at the offset location, the elastic interior is the
clear span, and **no new rigid-link code** is needed. `-offset offI offJ` is the
element-level convenience; under the hood it is the transformation's `-jntOffset`
along the member axis. Steel (offset to RBS centroid) and concrete (offset to
column face) differ only by the offset value.

---

## 4. Material interface — agnostic, tags only

The element takes **pre-built `uniaxialMaterial` moment-rotation tags** and is
indifferent to the law:

| use | uniaxial law |
|---|---|
| steel beam | `IMKBilin` / `Bilin` (tag 87 family) / `ModIMKPeakOriented` |
| concrete beam | `ModIMKPeakOriented` / pinching (`Pinching4`, `ModIMKPinching`) |
| any future law | whatever moment-rotation `uniaxialMaterial` you author |

"Both steel and concrete" is satisfied by **tag choice**, with zero element-side
branching. Raw IMK backbone parameters live in the uniaxial material's own (easy)
API, not in the element command.

### 4.1 Command

```tcl
element LadrunoIMKBeam $tag $iNode $jNode $transfTag \
    -E $E -G $G -A $A -J $J -Iy $Iy -Iz $Iz \
    -hinge both \                  ;# both | i | j
    -matZ  $imkZ_i $imkZ_j \       ;# strong-axis IMK tags (1 value = symmetric)
    -matY  $imkY_i $imkY_j \       ;# weak-axis (optional; omit → My elastic)
    -offset $offI $offJ            ;# column-face offsets (member-axis)
    [-shear $GAsy $GAsz]           ;# v1.x: Timoshenko interior (default off)
```

Python mirror via `ops.element('LadrunoIMKBeam', tag, i, j, transfTag, '-E', E, ...)`.

---

## 5. Recorders (single tag — the usability win)

| response | content |
|---|---|
| `hingeRotation` | `θ_h` per hinged end/axis (the IMK strains) |
| `hingeMoment` | `M` per hinged end/axis (the IMK stresses) |
| `basicForce` | `q = [N, Mz_i, Mz_j, My_i, My_j, T]` |
| `basicDeformation` | `v` |
| `material $end $axis ...` | passthrough to the IMK material (peak rotation, cyclic deterioration state, etc.) |

---

## 6. Validation plan

Zone-A pytest battery `tests/test_ladrunoIMKBeam_element.py`:

1. **Elastic degenerate mode** — no materials ⇒ matches `elasticBeamColumn`
   stiffness/eigenpairs to ~1e-9 (proves the basic assembly + transform).
2. **Single-hinge monotonic** — strong-axis cantilever, one IMK end; element
   moment-rotation reproduces the bare `uniaxialMaterial` backbone (the series
   with the elastic interior is the only addition; n-factor-free exactness).
3. **Lignos-Krawinkler cyclic** — steel MRF beam under a standard cyclic protocol
   vs a published `ModIMKPeakOriented`/`Bilin` backbone (cyclic deterioration,
   energy).
4. **Diaphragm sanity** — beam with a `rigidDiaphragm`-imposed spurious axial
   `P`; assert the moment-rotation law is bit-unchanged vs the no-`P` run (the
   core motivating guarantee) **and** the solver converges where the fiber-hinge
   `forceBeamColumn` does not.
5. **Column-face offset** — moment at the hinge vs moment at the joint differ by
   `V·offset`; assert against closed form.
6. **Python series-spring oracle** — the two-element (elastic + 2 `zeroLength`,
   n-factor) builder, used **only as a cross-check**, matches this element's
   global response to ~1e-6 on a fixed pushover. (Python is the oracle; the C++
   element is the deliverable.)

---

## 7. Authoring checklist (per `CLAUDE.md`)

- [ ] classTag **33003** in `SRC/classTags.h` (Element band) + this row in
      `LEDGER_implementations.md`.
- [ ] `FEM_ObjectBroker` registration (`getNewElement`).
- [ ] `OPS_LadrunoIMKBeam.cpp` parser; Tcl + Python wiring.
- [ ] LADRUNO header stamp — add the new files to `Ladruno_scripts/stamp_headers.py`
      GLOBS and rerun (non-optional for fork-authored classes).
- [ ] Banner line in `Ladruno_scripts/banner_features.txt` → `patch_banner.py`.
- [ ] Zone-A battery (§6); manifest/CI gates.
- [ ] `sendSelf`/`recvSelf` (rebuild basic assembly on receive; serialize the 4
      committed `θ_h` + material tags + elastic props + offsets + hinge mode).
- [ ] PR based on `ladruno`.

---

## 8. Open / deferred

- **Timoshenko interior** (`-shear`) — v1.x flag; default Euler-Bernoulli.
- **Biaxial-coupled column hinge** — `LadrunoIMKColumn` / `zeroLengthND` +
  coupled moment-surface `nDMaterial`; the explicit v1 non-goal (§1.2).
- **Axial/shear hinges** — the basic system has the slots (`N`, and shear via the
  transform); not exposed in v1.
- **End-moment release (true pin, `M ≡ 0`)** — *distinct from "elastic."* An
  **elastic** end (no material in that slot — §2.3) still carries moment with the
  full elastic rotational stiffness; a **released** end carries **zero** moment
  and frees the rotation, reducing the active end stiffness from `4EI/L` to
  `3EI/L` (the far end statically condensed out). There is no `-release` flag yet;
  it is a clean future addition (per-end-per-axis, exact condensation — a
  released end is "the hinge whose moment is identically zero").

  **Workaround available today (no code):** give that end a near-zero-stiffness
  `Elastic` hinge material — the series flexibility `1/k → ∞` drives the end
  moment to ~0. Use roughly `k_release ≈ 1e-5 · (4EI/L)` of the bending axis
  (i.e. ~1e-5 of the member's elastic rotational stiffness):

  ```python
  # strong-axis pin at end j (Iz axis): elastic interior + a ~zero-stiffness hinge
  ops.uniaxialMaterial('Elastic', 7, 1.0e-5 * (4.0*E*Iz/L))
  ops.element('LadrunoIMKBeam', tag, i, j, A,E,G,Jx,Iy,Iz, transf, '-matZj', 7)
  ```

  Verified: this gives `k_i = 3EI/L` to ~5 figures and residual `M_j/M_i ≈ 7e-6`
  — effectively a perfect pin. Stay above the element's `ktFloor` guard
  (`1e-8·4EI/L`); a factor in `[1e-6, 1e-4]` is the sweet spot (smaller hurts
  conditioning, larger leaves a non-negligible residual moment). See
  `LEDGER_quirks.md`.
- **Sensitivity** (`DDM`) — standard element plumbing, deferred.

Related: [[10_ladruno_j2_plasticity]], [[13_ladruno_uniaxial_j2_adr]] (candidate
fiber laws if a fiber-section hinge variant is ever wanted), and the
ZeroLength/SectionAggregator reference in
`Ladruno_implementation/zerolength_and_link_springs_guide.md` (the vanilla path
this element packages).
