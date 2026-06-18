---
title: LadrunoCohesiveHingeBiaxial — coupled Mz–My cohesive interaction surface
project: Ladruno
status: PR-4a shipped (coupled material + element wiring; B-K mode-mix / consistent off-radial tangent / torsion deferred)
priority: medium
owner: nmora
tags:
  - adr
  - material
  - cohesive
  - biaxial
  - embedded-discontinuity
  - frame
---

# ADR 34 — LadrunoCohesiveHingeBiaxial (coupled biaxial cohesive)

> Builds on [[33_ladruno_dispbeamcolumn3d_hinge_adr]] (PR-3b shipped the **biaxial
> element** carrying two **block-diagonal** scalar `LadrunoCohesiveHinge` laws; the
> dominant biaxial interaction already flows through the bulk section tangent
> `ks(MZ,MY)` in `K_αα`). This is the **v2 cohesive refinement** reserved there:
> a single coupled `M([[θ]])` law with an **Mz–My interaction surface**, so the
> hinge respects a biaxial moment capacity (an ellipse) instead of two independent
> axes. Generalizes the scalar [[LadrunoCohesiveHinge]] (`MAT_TAG 33003`).

## What

A new fork-authored **`NDMaterial` of order 2** — strain = the rotation-jump vector
`α = [α_z, α_y]`, stress = the cohesive moment vector `M = [M_z, M_y]`, tangent = the
2×2 `dM/dα`. A discrete cohesive law carrying **fracture energy per axis** (`Gf_z`,
`Gf_y`) directly (no `lch`), with an **elliptical activation surface** and **isotropic
secant damage** — a mixed-mode generalization of the scalar rigid-softening hinge.

Wired into `LadrunoDispBeamColumn3d` as a **third hinge mode** via `-hingeBiaxial
$ndMatTag` (mutually exclusive with `-hinge`/`-hingeY` and `-nl`). The element's
existing biaxial condensation (PR-3b: unified inner Newton + eigenvalue-floored 2×2
`K_αα` + rank-2 condensation) is reused **verbatim**; the only element change is that
the cohesive contribution to `K_αα` becomes a **full 2×2** (the coupled law's off-diagonal
`dM_z/dα_y`) instead of a diagonal `[Kz, Ky]`.

`ND_TAG_LadrunoCohesiveHingeBiaxial = 33004` (ND registry — free; the `33004` slot ADR 33
reserved in the uniaxial/nD band, here realized in the ND registry since the law is vector-valued).

## Formulation (pinned — reduces EXACTLY to the scalar law on each pure axis)

Per-axis penalty + guarded floor, identical to the scalar law applied per axis:
`Kpen_i = penaltyRatio · Mc_i²/(2 Gf_i)` (default ratio 1000), floored at `Mc_i²/(2 Gf_i)`;
activation jump `α0_i = Mc_i/Kpen_i`.

**Normalized jump / elliptical norm:** `e_i = α_i/α0_i`, `r = √(e_z² + e_y²)`.
- `r ≤ 1`: **closed** (rigid), `M_i = Kpen_i α_i`, no damage. At `r = 1` the moment lies
  on the activation ellipse `√((M_z/Mc_z)² + (M_y/Mc_y)²) = 1`.
- `r > 1`: damage onset. History `r_max = max r` (monotone, irreversible).

**Mode mix** (direction in normalized space): `w_z = e_z²/r²`, `w_y = e_y²/r²` (`w_z+w_y=1`).
Mixed peak-scale and fracture energy interpolate linearly (→ exact on each pure axis):
`S = w_z·Mc_z²/Kpen_z + w_y·Mc_y²/Kpen_y`, `Gf_mix = w_z Gf_z + w_y Gf_y`.
`S/2 < Gf_mix` always (per-axis floor), so a softening branch exists for every mix.

**Effective 1D law in `r`** (work-conjugate traction `T_eff`, `T_eff·dr = M·dα`):
`T_eff = (1−D)·S·r`; pre-peak `T_eff = S r` (D=0), softening calibrated to `∫T_eff dr = Gf_mix`:
- `EXP`: `T_eff = S·exp(−A(r−1))`, `A = S/(Gf_mix − S/2)`.
- `LINEAR`: `T_eff = S(1 − (r−1)/r_f)`, `r_f = 2(Gf_mix − S/2)/S`.

**Isotropic secant damage** `D = 1 − T_eff(r_max)/(S r_max)` (with `S, Gf_mix` at the current mix).
**Stress** `M_i = (1−D)·Kpen_i·α_i` (secant to origin → no strength recovery, no residual at α=0).

**Reduction check (pure-z, w_z=1):** `S = Mc_z²/Kpen_z`, `Gf_mix = Gf_z`, `M_z = T_eff·Kpen_z/Mc_z`
→ peak `Mc_z` at `r=1`, EXP rate `A = aExp_z·α0_z` (the scalar `aExp_z` mapped to normalized `r`).
**Identical** to `LadrunoCohesiveHinge(Mc_z, Gf_z)`. Same for pure-y and LINEAR. ∎

**2×2 tangent.** Loading (`r = r_max` increasing): `∂M_i/∂α_j = (1−D)Kpen_i δ_ij − Kpen_i α_i ∂D/∂α_j`,
with `∂D/∂α_j` taken through `r` at frozen mix (`∂r/∂α_j = α_j/(α0_j² r)`) — **exact on radial paths**
(the gates), a consistent mixed-mode tangent off-radial. Unloading (`r < r_max`): `D` frozen,
`∂M_i/∂α_j = (1−D)Kpen_i δ_ij` (diagonal secant). The `−Kpen_i α_i ∂D/∂α_j` term is the **cohesive
off-diagonal** that, added to the bulk `ks(MZ,MY)` term, makes the element's `K_αα` the true coupled 2×2.

## Where

- **New material:** `SRC/material/nD/LadrunoCohesiveHingeBiaxial.{h,cpp}` (+ `CMakeLists.txt`),
  registered at the 5 nD sites: `classTags.h`, `FEM_ObjectBrokerAllClasses.cpp`
  (`getNewNDMaterial` case), `OpenSeesNDMaterialCommands.cpp` (Python), the Tcl
  `nDMaterial` table, `CMakeLists.txt`.
- **Element wiring:** `SRC/element/ladrunoDispBeamColumn/LadrunoDispBeamColumn3d.{h,cpp}` +
  `OPS_LadrunoDispBeamColumn.cpp` — `-hingeBiaxial $matTag` holds an `NDMaterial *theHingeC`;
  a coupled branch of the biaxial inner Newton uses the full 2×2 cohesive tangent. Serialization
  grows for the nD-material send (broker `getNewNDMaterial`). Gated so `-hinge`/`-hingeY` and the
  no-hinge path are byte-identical.
- **Tests:** `tests/test_ladrunoCohesiveHingeBiaxial_material.py` (material-point: on-axis reduction
  to the scalar law, elliptical onset, mixed-mode dissipation, irreversibility, tangent FD, roundtrip)
  + biaxial-element gates in `tests/test_ladrunoDispBeamColumn3d_hinge_biaxial.py` (coupled hinge
  dissipates on a combined path; reduces to the block-diagonal element on pure axes).

## Decisions / scope

- **Linear mode-mix** `Gf_mix = w_z Gf_z + w_y Gf_y` (B-K/power-law variants deferred); radial-path
  exact, instantaneous-mix off-radial (documented v1 approximation, like most engineering CZMs).
- **NDMaterial order 2**, not a bespoke class — reusable, broker-serializable, standard plumbing.
- **No torsion coupling** (the twist channel stays linear; `α_t` is ADR 33 PR-3c).
- The element keeps `-hinge`/`-hingeY` (block-diagonal) as the cheaper default; `-hingeBiaxial`
  is the opt-in coupled law for cases needing the moment interaction ellipse.

## Implementation log

### 2026-06-18 — PR-4a SHIPPED: LadrunoCohesiveHingeBiaxial + element wiring

The coupled law and its element wiring, implemented and validated.

- **Files:** `SRC/material/nD/LadrunoCohesiveHingeBiaxial.{h,cpp}` (`ND_TAG 33004`) +
  registration at the 4 nD sites (`classTags.h`, `FEM_ObjectBrokerAllClasses.cpp`,
  `OpenSeesNDMaterialCommands.cpp`, `material/nD/CMakeLists.txt`);
  `SRC/element/ladrunoDispBeamColumn/LadrunoDispBeamColumn3d.{h,cpp}` +
  `OPS_LadrunoDispBeamColumn.cpp`; `tests/test_ladrunoCohesiveHingeBiaxial`… (in
  `tests/test_ladrunoDispBeamColumn3d_hinge_coupled.py`, 9 tests).
- **Material:** order-2 NDMaterial, exactly the pinned formulation (elliptical onset, mixed-mode
  isotropic secant damage, per-axis Gf). `getEnergy()` = 2D path work; scalar responses
  (energy/momentZ/momentY/jumpZ/jumpY/rMax/damage). Reduces to `LadrunoCohesiveHinge` on each axis.
- **Element:** `-hingeBiaxial $ndMatTag` holds `NDMaterial *theHingeC` (exclusive with
  `-hinge`/`-hingeY`/`-nl`). The PR-3b biaxial inner Newton (`solveHingeJumpBiaxial`) is REUSED — it
  branches the cohesive evaluation: one coupled law returning a moment vector + a **full 2×2**
  cohesive tangent (vs the two scalar laws' diagonal), so the cohesive off-diagonal `Czy` enters
  `K_αα`. Same eigenvalue-floored inverse + rank-2 condensation. Serialization grows `data` 23→26
  (the nD material's classTag/dbTag ride in the data Vector — folded in to avoid an FE_Datastore
  ID-size collision — and it sendSelf's via `getNewNDMaterial`). `"hingeBiaxial <resp>"` passthrough.
- **Inner-Newton robustness (the one correction):** the coupled tangent is sign-discontinuous at the
  elliptical onset (r=1) and frozen-mix (exact only on radial jump paths), so a full Newton step can
  overshoot the activation kink. Added **adaptive step damping** (halve when the residual grows,
  relax back otherwise) + `maxIter` 30→50. It is a NO-OP while the residual decreases monotonically,
  so the block-diagonal PR-3b path is unchanged (28+15 element/biaxial tests still bit-stable).
- **Verified (OpenSeesPy, 9/9; full 88/88 element+hinge+cohesive 2D+3D):** parse + exclusivity;
  **on-axis reduction** — pure-Mz coupled element matches the scalar `-hinge` element's Mz path to
  1e-6 and dissipates `Gf_z` with the weak axis closed (symmetric pure-My → `Gf_y`); **elliptical
  interaction onset** — a proportional prescribed path activates on `√((Mz/Mcz)²+(My/Mcy)²)=1`,
  strictly inside the independent capacities; **mixed-mode dissipation** between `Gf_y` and `Gf_z`;
  **coupled tangent** — product-of-inertia section + coupled hinge through deep biaxial softening
  converges under tight Newton; **finite-rotation invariance** under `CorotCrdTransf3d` (>0.5 rad,
  both jumps open); DB roundtrip both-open.
- **Honest scope (v1):** deeply weak-axis-dominant paths driven by INDIRECT control (a non-radial
  jump path that lands the global solver on the unstable elliptical onset) sit at the edge of the
  frozen-mix tangent; prescribe the rotations (radial control) there. The cohesive mode-mix is the
  instantaneous direction (radial-exact). **NEXT:** B-K/power-law mode-mix; a fully consistent
  off-radial tangent (∂mix/∂α); torsion (ADR 33 PR-3c).
