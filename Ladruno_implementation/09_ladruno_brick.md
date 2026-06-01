---
title: LadrunoBrick — unified hexahedral solid element
project: Ladruno
status: draft
priority: high
owner: nmora
tags:
  - implementation
  - element
  - solid
  - hexahedron
  - locking
---

# LadrunoBrick — unified hexahedral solid element

## What

A single 8-node hexahedral solid element that exposes its **formulation** as a
parameter instead of forcing the user to pick a different element class. One
class, one `classTag` (**`ELE_TAG_LadrunoBrick = 33002`**, next in the Ladruno
element band after `BezierTri6=33000`/`BezierTet10=33001`), one Tcl/Python
command, with the anti-locking treatment chosen at construction:

```
element('LadrunoBrick', tag, n1..n8, matTag, '-formulation', <std|bbar|uri|eas>, ...)
```

This collapses the upstream `Brick` / `bbarBrick` split (and the *missing*
reduced-integration explicit hex, and the *missing* enhanced-strain hex) into
one verifiable family.

**Scope — v1 (this plan):** small-strain (geometrically linear) kinematics, and
the first three formulations:

- `std`  — full 2×2×2 Gauss, displacement (reproduces upstream `Brick`)
- `bbar` — mean-dilatation B-bar (reproduces upstream `bbarBrick`)
- `uri`  — 1-point reduced integration **+ hourglass control** (the cheap
  explicit workhorse; new to OpenSees' brick family)

**Scope — explicitly deferred:**

- `eas` (enhanced assumed strain) → **v2**. It carries internal element state
  `α` (static condensation, commit/`sendSelf` of `α`) that dwarfs the other
  three combined; land the cheaper wins first. The `-formulation eas` keyword is
  *reserved* in v1 (errors with "not yet implemented") so the API is stable.
- Large displacement / large rotation (`corot`) and large strain (`finite`) →
  handled by a separate, orthogonal **geometry-method layer**, planned in
  [[solid_transformation_wrapper]]. v1 is the `linear` geometry only, but its
  kernel is carved so that layer drops in without a rewrite (see *Three seams*).

> **Design rule.** `-formulation` is a *single selector*, not four independent
> booleans. Most boolean combinations are redundant or contradictory:
> **B-bar ≡ selective ("partial") integration of the volumetric term** (Hughes
> 1980 — the same element stiffness, two derivations), and `bbar` presupposes a
> *fully* integrated deviatoric part so it cannot combine with `uri`. A selector
> makes invalid states unrepresentable.

## Why

The upstream 8-node `Brick` ([SRC/element/brick/Brick.cpp](../SRC/element/brick/Brick.cpp))
is a correct textbook trilinear hex, but it carries the textbook *limitations*:

1. **Volumetric locking** near-incompressible (ν→0.5, J2 at the isochoric yield
   surface, saturated/undrained soil) — full integration over-constrains
   `∇·u=0`. The cure (`bbar`) lives in a *separate class* today, so a model must
   be rebuilt to switch.
2. **Shear / parasitic-bending locking** under bending-dominated regimes — *not*
   cured by B-bar. There is no enhanced-strain or reduced-integration hex in the
   brick family to fix it.
3. **No cheap explicit hex.** Only full 2×2×2 exists. Explicit dynamics — the
   fork's direction ([[05_robust_central_difference]], [[Ladruno_explicit_roadmap]])
   — wants a 1-point hex with hourglass control (~8× cheaper internal-force
   evaluation, the LS-DYNA workhorse).
4. **A real bug.** `Brick::setParameter` loops `for (i=0;i<4;i++)` over
   `theDamping[i]` but the array has **8** entries
   ([Brick.cpp:1962](../SRC/element/brick/Brick.cpp)) — damping objects 4–7 never
   receive the parameter (copy-paste from the 4-node quad). We won't inherit it.

We want one solid element that is **verifiable** (reproduces the upstream
elements bit-for-bit where they overlap), **explicit-ready** (non-negative
lumped mass, reduced-integration option), **recorder-native** (honours the
[[ladruno_element_contract]] so the Ladruno recorder captures it with zero
recorder edits), and **future-proof** (large-rotation/large-strain drop in as a
geometry layer, not a fork).

## Where

- **New code:**
  - `SRC/element/ladrunoBrick/LadrunoBrick.{cpp,h}` — the element.
  - `SRC/element/ladrunoBrick/TclLadrunoBrickCommand.cpp` — Tcl command.
  - `SRC/element/ladrunoBrick/CMakeLists.txt`, `Makefile`.
  - `tests/test_ladrunoBrick_element.py` — Zone-A battery (see *Testing*).
- **Modify (vanilla — log in [[LEDGER_vanilla_files]] with `// Ladruno` tags):**
  - `SRC/classTags.h` — `#define ELE_TAG_LadrunoBrick 33002`.
  - `SRC/domain/component/FEM_ObjectBroker.cpp` — broker entry for `recvSelf`.
  - `SRC/interpreter/OpenSeesElementCommands.cpp` (+ Tcl `TclElementCommands`) —
    register `OPS_LadrunoBrick`.
  - `Ladruno_scripts/banner_features.txt` — one line (then `patch_banner.py`).
- **Reference (copy patterns):**
  - `SRC/element/brick/Brick.cpp` — baseline state cycle, mass, body force,
    `setResponse` stress/strain tree, `sendSelf`/`recvSelf`.
  - `SRC/element/brick/BbarBrick.cpp:1098` — the dev/vol B-bar split.
  - `SRC/element/brick/shp3d.cpp` — trilinear shape functions + Jacobian.
  - `SRC/element/bezierTetrahedron/BezierTet10.cpp` — sibling Ladruno solid:
    non-negative lumped mass, `|detJ|` orientation-robustness, recorder probes.

## How

### Public API

```
element('LadrunoBrick', eleTag, n1,n2,n3,n4,n5,n6,n7,n8, matTag
        [, '-formulation', <std|bbar|uri|eas>]   # default: std
        [, '-hourglass', <type>, <coeff>]        # uri only; type ∈ {viscous,stiffness,physical}
        [, '-lumped']                            # lumped mass (explicit); default consistent
        [, '-b', bx, by, bz]                     # body force
        [, '-damp', dampTag])
```

Tcl mirror: `element LadrunoBrick $tag $n1 ... $n8 $matTag -formulation bbar ...`.

### Internal structure — the formulation selector

```cpp
enum class Formulation { STD, BBAR, URI, EAS };   // EAS reserved (v2)
enum class GeomMethod  { LINEAR, COROT, FINITE };  // v1: LINEAR only (see seam 2)

// dispatch in formResidAndTangent():
//   STD  -> B from computeBLinear  (8-pt Gauss)
//   BBAR -> B from computeBbar     (8-pt dev + 1-pt mean dilatation)
//   URI  -> B at centroid (1-pt) + addHourglassStabilization()
//   EAS  -> error in v1
```

`std`/`bbar` share Brick/BbarBrick's exact kernels (we *copy*, not subclass, to
own the bug-fix and the seams). `uri` is the new path.

### Reduced integration + hourglass control (`uri`)

1-point Gauss captures only the 6 constant-strain modes + 6 rigid-body modes of
the 24-DOF element → **12 hourglass (zero-energy) modes** remain and must be
stabilized. Use the Flanagan–Belytschko (1981) hourglass base vectors `Γ_α`
(α=1..4), 3 directions → 12 modes. Two stabilization flavours, picked by target
solver (we target **both** implicit and explicit):

- **`viscous`** / **`stiffness`** (Flanagan–Belytschko) — a user coefficient
  scales a perturbation hourglass force/stiffness. Cheap; standard for
  **explicit**. `viscous` damps hourglass velocity; `stiffness` adds an elastic
  hourglass restoring force.
- **`physical`** (Belytschko–Bindeman assumed-strain / Puso) — consistent
  assumed-strain stabilization sampled from the centroid material tangent; **no
  tuning knob**, consistent tangent for **implicit/static**. Default for
  implicit.

Default: `physical` when available for the material's tangent, else
`stiffness` with coeff 0.05. Document the knob; never silently zero it.

### Mass

Reuse Brick's consistent/lumped paths. For explicit, lumped must be
**non-negative** per Gauss/node (trilinear hex row-sum lumping is already
non-negative — unlike higher-order/serendipity; cf. [[06_bezier_tet10]]'s
`ρVe/10` care). Mass is **formulation-independent** (always full integration of
∫ρNᵀN) — do *not* lump-integrate mass with the URI rule.

### The three seams (carve in v1, even though v1 only uses the first state)

This is the load-bearing design decision — it makes the "beast" roadmap additive
instead of a rewrite. See [[solid_transformation_wrapper]] for seams 2–3 in full.

| Seam | What it isolates | v1 behaviour | Later |
|---|---|---|---|
| **1. Kinematics ledger** | compute strain measure from (X, u, ∇N) at a GP — *intrinsic, present in every method* | emits engineering ε = B·u | also emit **F** = I + Σ uₐ⊗∇Nₐ |
| **2. Geometry method** | de-rotate-in / re-rotate-out around an oblivious kernel | identity (`linear`) | `corot` (polar R + K_geo), `finite` (TL) — [[solid_transformation_wrapper]] |
| **3. Material adaptor** | element↔material boundary | pass ε straight to `setTrialStrain` | log-strain adaptor: F→E_log→return-map→τ,ℂ (**owned by the parallel material-wrapper work**) |

The constitutive/integration **core** (the `bbar`/`uri`/`eas` logic) sits
*inside* seams — it never learns which geometry method or material adaptor wraps
it. That orthogonality is why a future `LadrunoBrick -formulation uri -geom corot`
needs no new kernel code.

### Recorder contract

A standard Lagrange hex is covered by the recorder's legacy class-tag fallback,
so `basisInfo` is **optional** ([[ladruno_element_contract]] §0 precedence). v1
implements the Part C result tree (`stresses`/`strains` per-GP, matching Brick's
`setResponse`) and the `stress3D6`/`strain3D6` averaged responses for vtkhdf.
When seam 2 gains `corot`, set `frameTimeVarying=1` and answer `localAxes` with
the current frame (contract §B.2).

### Testing — Zone-A battery (`tests/test_ladrunoBrick_element.py`)

The verification ladder, gated in CI (Zone-A, numpy-based like the CDL/Bathe
batteries):

1. **Regression (the headline gate):** `-formulation std` reproduces upstream
   `Brick` tangent & resisting force to ~1e-12 on a distorted hex; `-formulation
   bbar` reproduces `bbarBrick`. (Bit-for-bit overlap is our correctness anchor —
   only possible because v1 is small-strain.)
2. **T0 patch test:** constant-strain patch on a distorted mesh, all formulations.
3. **Rank / spectral:** stiffness rank = 24−6 = 18 for `std`/`bbar`; for `uri`
   confirm 12 hourglass modes are stabilized back to rank 18 with **no spurious
   energy** in the hourglass patterns.
4. **Volumetric locking:** near-incompressible (ν=0.4999) thick-walled / confined
   block — `std` locks (tip displacement ≪ analytic), `bbar`/`uri` recover it.
5. **Bending / shear locking:** slender cantilever in pure bending — `std`/`bbar`
   over-stiff, `uri` (and later `eas`) near-exact.
6. **Explicit:** `dt_cr` sanity + a wave-propagation / Cook's-membrane run under
   `CentralDifferenceLadruno` with `-formulation uri -lumped` (no hourglass
   growth).
7. **Recorder round-trip:** record `stresses`/`strains` via the Ladruno recorder,
   read back, assert shape/values.

## Risks / open questions

> [!question]
> **Hourglass default for `uri`.** `physical` (assumed-strain) is the right
> default but needs the centroid material tangent — fine for elastic/J2, but does
> a generic `NDMaterial::getTangent` at the centroid give a usable stabilization
> modulus for all materials? Fallback to `stiffness`+coeff if not. Decide per
> material class or globally?

> [!question]
> **`std`/`bbar`: copy vs. share with upstream Brick.** Copying lets us fix the
> `setParameter` damping bug and carve the seams, at the cost of duplicated
> kernel. Subclassing `Brick` would inherit the bug and the seam-less structure.
> Leaning **copy** (own the code) — confirm.

- **EAS deferral discipline:** keep `-formulation eas` parsing but error clearly;
  do not half-build the `α` state path in v1 (it touches `commitState`,
  `revertToLastCommit`, `sendSelf`/`recvSelf`).
- **Material adaptor coordination:** seam 3's signature must match what the
  parallel material-wrapper work expects (F in, S/ℂ or τ/c out). Lock the
  boundary in [[solid_transformation_wrapper]] before either side hardens.
- **Lumped mass for explicit:** trilinear hex is safe; if v-next adds 20/27-node
  variants, revisit non-negativity (HRZ vs row-sum).
- **Backwards compat:** new class/tag — no impact on existing models. The upstream
  `Brick`/`bbarBrick` stay; `LadrunoBrick` is additive.

## Roadmap

- **v1 (this plan):** small-strain `std`/`bbar`/`uri`(+hourglass), 8-node, full
  test ladder, recorder-native, banner + ledgers. Bug-free `setParameter`.
- **v2:** `eas` (enhanced assumed strain, internal `α` + static condensation);
  `corot` geometry method via [[solid_transformation_wrapper]] (large rotation,
  reuses every v1 formulation + the full small-strain material library).
  **`eas` blueprint = OpenSees `SSPbrick`** (`SRC/element/UWelements/SSPbrick.cpp`):
  it is exactly `bbar` (constant part via mean-dilatation `dNmod`, `:1266`) **+
  statically-condensed EAS** (9 enhanced modes, `interior = FCF − K_uα K_αα⁻¹ K_αu`,
  `:1968`). That condensation is what makes it accurate for **all ν** — a fixed
  assumed strain (our `physical`) cannot. Cross-validate v2 against `SSPbrick`.
- **v3 (if geotech large-deformation / base-isolation drives it):** `finite`
  (total-Lagrangian) geometry + the log-strain **material adaptor**, which
  promotes the existing small-strain soil/metal materials to large strain.
  Possibly higher-order (20N serendipity / 27N Lagrange) siblings.

## Implementation log

- 2026-05-31 — Plan drafted. Single `-formulation {std|bbar|uri|eas}` selector
  (B-bar ≡ selective/partial integration → one name, not two flags; bbar⊥uri).
  v1 = small-strain std/bbar/uri+hourglass; EAS→v2; corot/finite→
  [[solid_transformation_wrapper]] geometry layer. Three seams carved
  (kinematics ledger / geometry method / material adaptor) so large-rot &
  large-strain are additive. classTag 33002 reserved. Headline gate: reproduce
  upstream `Brick`/`bbarBrick` to ~1e-12. Noted the upstream `setParameter`
  damping `i<4`-over-8 bug we will not inherit.
- 2026-06-01 — **SHIPPED v1: std + bbar + uri (stiffness + physical hourglass).**
  std↔`Brick` / bbar↔`bbarBrick` to ~1e-9 (PR #69, merged). uri = 1-pt + FB
  `stiffness` hourglass (PR #69). `physical` = full-normal assumed strain +
  Belytschko 8.7.26 mode-subset shear, validated by a bending-convergence
  benchmark cross-checked vs `SSPbrick` (PR #72). **physical is a SHEAR-locking
  cure only** (matches SSPbrick at ν=0, converges 0.94→1.005) — it **volumetric-
  locks at ν→0.5** (use `bbar`). Finding: the eq-8.7.26 isochoric dev-projection
  over-softens bending (worse with ν); no fixed assumed strain cures both
  shear+volumetric across ν. The general-ν element needs **EAS** (`eas` v2) —
  and `SSPbrick` is the blueprint (bbar + condensed EAS, `:1968`). See
  [[LEDGER_quirks]]. Composes with the `SolidTransformation` seam (#71): physical
  routes its strain through `computeLocalDisp()` (uCore). Reserved still: `eas`,
  `uri`+`viscous`. Banner line still pending (ships when merged).
