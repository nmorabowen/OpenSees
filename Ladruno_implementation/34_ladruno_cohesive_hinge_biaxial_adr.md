---
title: LadrunoCohesiveHingeBiaxial — coupled Mz–My cohesive interaction surface
project: Ladruno
status: PR-4b shipped (coupled material + element wiring + B-K mode-mix + consistent off-radial tangent; torsion deferred)
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

- **Mode-mix** `Gf_mix` = Benzeggagh-Kenane `Gf_z + (Gf_y - Gf_z) w_y^eta` via `-bk eta` (default
  `eta=1` = the original linear interpolation, bit-identical). The **consistent off-radial tangent**
  (PR-4b) carries `dD/dw . dw/dalpha`; see the PR-4b log for the mix-invariance subtlety that makes
  it a no-op at `eta=1`.
- **NDMaterial order 2**, not a bespoke class — reusable, broker-serializable, standard plumbing.
- **No torsion coupling** (the twist channel stays linear; `α_t` is ADR 33 PR-3c).
- The element keeps `-hinge`/`-hingeY` (block-diagonal) as the cheaper default; `-hingeBiaxial`
  is the opt-in coupled law for cases needing the moment interaction ellipse.

## Implementation log

### 2026-06-19 — PR-4b SHIPPED: B-K mode-mix + consistent off-radial tangent (item #1 + #2)

Closing handoff items #1 (consistent off-radial tangent) and #2 (B-K mode-mix) together — they
turned out to be inseparable.

**The load-bearing discovery (don't re-derive — it redefined item #1).** The v1 "frozen-mix tangent,
exact only on radial paths" caveat is **vacuous for the default linear mix**: the exact-reduction
calibration `Kpen_i = ratio*Mc_i^2/(2 Gf_i)` makes `c_i = Mc_i^2/Kpen_i = 2 Gf_i/ratio`, so
`S = (2/ratio) Gf_lin` and `Esoft = Gf_lin(1 - 1/ratio)` — both proportional to the LINEAR `Gf_lin`,
which makes the softening rate `A = S/Esoft = 2/(ratio-1)` a **constant, independent of the mode mix**.
Hence `D = D(r)` ALONE (verified against the compiled material: at fixed `r`, damage spread across
pure-z / 45deg / pure-y is `0` for EXP, `1e-16` for LINEAR). The v1 tangent was therefore ALREADY the
exact consistent off-radial Jacobian; the `dD/dw . dw/dalpha` term item #1 asked for is identically
zero at `eta=1`.

**So #1 only becomes real WITH #2.** A **Benzeggagh-Kenane** fracture-energy mix
`Gf_mix = Gf_z + (Gf_y - Gf_z) w_y^eta` (`-bk eta`, default 1) breaks the `S ∝ Gf_mix` proportionality
for `eta != 1`, so `A` and hence `D` become genuinely mix-dependent and the consistent tangent's
mode-mix term goes LIVE and nonzero. Peak-scale `S` and the elliptical onset stay LINEAR/unchanged →
exact pure-axis reduction (`Gf_mix = Gf_z` on the z-axis for ANY `eta`) is preserved.

- **Material:** added `bkEta` (parsed `-bk`/`-eta`/`-modemix`, `>0`); `effectiveLaw` uses the B-K
  `Gf_mix` (branched on `eta==1.0` so the default is bit-for-bit unchanged). `setTrialStrain` adds the
  consistent mode-mix tangent term `dD/dw_z . dw_z/dalpha_j` (guarded to `eta!=1` AND both axes active
  — `pow(w_y, eta-1)` would blow up at `w_y=0` for `eta<1`, but `dw_z/dalpha` is zero there anyway).
  `dD/dw_z = -(1/r) d(rho)/dA . dA/dw_z` with `rho = T_eff/S`; `dGf_mix/dw_z = eta(Gf_z-Gf_y) w_y^(eta-1)`.
  Serialization `data` 13→14 (carry `bkEta`); JSON/Print show it; new raw `tangent` material response
  (id 207, `[K00,K01,K10,K11]`).
- **Magnitude caveat:** at the **near-rigid default `ratio=1000`** the mode-mix term is real but tiny
  (~1e-4 of the tangent) because `D ≈ 1 - 1/r` (mix-independent secant geometry) dominates and the
  mix only enters via `A = O(1/ratio)`. It is numerically significant only for softer penalties
  (~1e-2 at `ratio=10`). The element default is unaffected; `-bk` is an opt-in refinement.
- **Tests:** `tests/test_ladrunoCohesiveHingeBiaxial_material.py` (46) — a `zeroLengthND` material
  point driven to exact prescribed `alpha`, reading the raw 2x2 via the new response. The capstone
  `test_consistent_tangent_matches_fd_offradial` (ratio=10, eta in {0.5,1,2}) asserts C++ tangent ==
  central-FD of `M(alpha)` to ~1e-7 AND that the **frozen-mix tangent deviates from FD for eta != 1**
  (the term is necessary). Plus mix-invariance-only-at-eta=1, pure-axis-unaffected-by-eta, B-K
  dissipation `Gf_z + (Gf_y-Gf_z)(1/2)^eta`, pre-peak diagonal penalty. Element-level (in
  `test_ladrunoDispBeamColumn3d_hinge_coupled.py`): a `-bk 2` coupled element converges through deep
  biaxial softening; a pure-axis `-bk 3` element matches the `eta=1` element's Mz path bit-for-bit.
- **Verified:** full battery **136/136** (the 88 stay bit-stable → `eta=1` default unchanged;
  +2 element B-K, +46 material-point); `check_manifest` / `check_classtags` OK (no new classTag).

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

### 2026-06-18 — Adversarial review (32-agent workflow) dispositions

A find→perspective-diverse-verify→synthesize workflow over PR-3b + PR-4a raised 8 surviving
findings. Dispositions after checking each against the source (88/88 stays green):

- **"Stray brace truncates the parser" (claimed must-fix) — FALSE POSITIVE.** Line 197 `}` is the
  necessary close of the option `while` loop; flags after `-hingeBiaxial` ARE parsed in later
  iterations (the `-hingeBiaxial $tag -nl` exclusivity test passes precisely because `-nl` is parsed
  and then rejected). The finding's own structural description contradicts its conclusion; the
  verifiers pattern-matched "double brace = bug". No change.
- **"Inner Newton returns success on maxIter" (claimed must-fix) — kept BEST-EFFORT by design.**
  Returning `-1` to cut the global step (a literal reading of the ADR-33 §invariant note) was tried
  and **reverted**: a transient inner miss at the coupled-law elliptical onset (r→1) happens for
  trial `v` the global Newton probes mid-iteration, and aborting the step there breaks softening
  robustness AND conflicts with the OpenSees internal-iteration convention (`ForceBeamColumn` et al.
  return best-effort + warn). The global `NormDispIncr` test is the real arbiter; the converged
  global state always carries a converged inner state. Comments now document this. Robustly meeting
  the strict "cut the step" intent needs the deferred onset-tangent hardening (a consistent
  off-radial / line-searched inner solve), tracked under NEXT.
- **Eigenvalue floor omits the off-diagonal `Czy` (latent) — left diagonal-only ON PURPOSE.** Folding
  `|Kaa_zy|` into the floor over-floors the modes near the elliptical onset and destabilizes the
  radial inner Newton (broke 4 tests when tried). The diagonal scale is the validated choice.
- **`hingeMscaleZ/Y` lifetime ratchet (should-fix) — not a defect.** The running scale anchors the
  inner tolerance to ~`Mc` (the hinge caps the moment at `Mc`, so the peak it tracks ≈ `Mc`), giving
  `tol ≈ 1e-11·Mc` always — tight, not "grossly loose". The reviewer's per-solve reset would
  REINTRODUCE the tol-collapse for a fully-broken hinge (`M→0`) that the member was added to fix.
- **Stale condensation cache after `revertToLastCommit` (should-fix) — documented v1 scope.**
  Real latent hazard only off the v1 contract (full Newton + DisplacementControl, where `update()`
  always precedes `getTangentStiff`); 2D PR-2b empirically confirmed the residual is always post-
  `update()`. A staleness guard is deferred (needs a re-solve in `getTangentStiff`).
- **Nits applied:** corrected the `Esoft>0` comment to credit the `penaltyRatio≥1000` guard (not the
  per-axis floors); documented the coupled tangent's by-design row asymmetry (`Kc01≠Kc10` when
  `Gfz≠Gfy` — the correct Jacobian of a non-potential damage model; the element symmetrizes before
  the 2×2 inverse).
- **Dismissed on the mechanics** (refuted by verification, not resurrected): dD/dr sign, off-diagonal
  coupling sign, trapezoidal work integral, r=0 division guard, eigenvalue-floor magnitude.
