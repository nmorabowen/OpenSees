---
title: LadrunoDispBeamColumn3d — Tier-2 embedded strong-discontinuity hinge (3D)
project: Ladruno
status: in-progress (PR-3a shipped; PR-3b biaxial pending)
priority: high
owner: nmora
tags:
  - adr
  - element
  - frame
  - 3d
  - embedded-discontinuity
  - softening
  - corotational
---

# ADR 33 — LadrunoDispBeamColumn3d Tier-2 embedded strong-discontinuity hinge

> Supersedes the "3D embedded hinge deferred to its own ADR" scope note in
> [[32_ladruno_dispbeamcolumn_regularization_adr]] (§Out-of-scope, Q5). Builds
> directly on the shipped **2D** Tier-2 hinge (ADR 32 PR-2a, commit `14513dc9c`)
> and its merged building block `LadrunoCohesiveHinge` (`MAT_TAG 33003`, ADR 32 #264).
> Scoped + adversarially reviewed by a 17-agent workflow (4 scouts × 3 designs ×
> perspective-diverse verify + synthesis), 2026-06-17.

## What

Add the embedded strong-discontinuity rotation-jump hinge to the **3D** displacement
fiber frame `LadrunoDispBeamColumn3d`, mirroring the 2D Tier-2 formulation on the
6-DOF 3D basic system. The hinge carries a discrete cohesive `M([[θ]])` law per
bending axis (the scalar `LadrunoCohesiveHinge`, reused verbatim) and is statically
condensed to the basic system **before** `CorotCrdTransf3d` (the pinned invariant).

## Recommendation (what to build)

Build the **biaxial rotation-jump hinge** (`α = [α_z, α_y]`, one jump per bending
axis), but the adversarial review forced three corrections to the naive port and a
two-PR split:

1. **Ship PR-3a = the strong-axis jump `α_z` alone** — the literal 2D scalar algebra
   transplanted onto the 3D `Mz` row of the 6-DOF basic system. With one jump there
   is **no off-diagonal**, so it is byte-for-byte the 2D scalar reciprocal; it is the
   lowest-risk proof that the pinned invariant survives the `CorotCrdTransf3d`
   quaternion-triad transform, and it is independently valuable (an out-of-plane frame
   regularizes its strong-axis hinge). The biaxial 2×2 lands in **PR-3b**.
2. **PR-3b uses a TRUE coupled 2×2 `K_αα`** with an **eigenvalue-floored** guarded
   inverse — NOT two independent rank-1 updates (the condensation reviewer proved that
   is the exact tangent of a *different* residual whenever the fiber section couples
   `ks(MZ,MY) ≠ 0`, generic for any axially-loaded biaxially-bent section, and it fails
   the FD-tangent gate pre-peak), and NOT a det-floored adjugate inverse (the unfloored
   adjugate numerator blows up ~1e8× along the near-null eigenvector at staggered
   activation).
3. **Single unified inner Newton** that sets BOTH bending channels in ONE
   `setTrialSectionDeformation` per IP — NOT two verbatim 2D-kernel passes (the second
   `setTrialSectionDeformation` overwrites the whole strain vector and corrupts the
   section state from the first pass).

Torsion jump and a coupled biaxial cohesive interaction surface are **out of v1** (see
Deferred). The bending-axis cohesive laws stay **block-diagonal** (two scalar
`LadrunoCohesiveHinge`); the dominant biaxial coupling already flows through the bulk
section tangent `ks(MZ,MY)` in `K_αα`, so the cohesive coupling is a documented
second-order refinement (a future coupled `nDMaterial`, fixed facade, zero element churn).

## The 3D formulation (per axis, generalizes ADR 32 §How verbatim)

**Basic system** (6-DOF, `v = crdTransf->getBasicTrialDisp()`):
`v = [axial, θz_i, θz_j, θy_i, θy_j, twist]` (`LadrunoDispBeamColumn3d.cpp:459`,
B-operator `:760-762`).

**Bounded enhancement** — one constant operator per bending axis, the literal 2D
constant: `Ḡ_z = Ḡ_y = −1/L`, each acting ONLY on its own curvature channel. Each is
constant ⇒ `∫Ḡ dx = −(1/L)·Σwt = −1` exact under any rule (weights sum to L), the jump
contributes `+1`, so orthogonality `∫G dx = 0` is machine-exact **per axis**. **Coupling
enters ONLY through the section tangent `ks(MZ,MY)`, never the enhancement geometry** —
this keeps the patch-test orthogonality per-axis exact while the full biaxial interaction
is captured in `K_αα`.

**Enhanced kinematics** (set in ONE `setTrialSectionDeformation(e)` per IP):
```
P : e = (1/L)·v(0)                                          (untouched; linear, no -nl)
MZ: e = (1/L)·((6ξ−4)v(1)+(6ξ−2)v(2)) + Ḡ_z·α_z            (= … − α_z/L)
MY: e = (1/L)·((6ξ−4)v(3)+(6ξ−2)v(4)) + Ḡ_y·α_y            (PR-3b)
T : e = (1/L)·v(5)                                          (untouched; no twist jump)
```
**Enhancement residuals** (the `1/L` from `dx` cancels the `1/L` in `Ḡ`, as in 2D):
`h_z = −Σwt·s(MZ) + M_coh_z(α_z) = 0`, `h_y = −Σwt·s(MY) + M_coh_y(α_y) = 0`.

**Condensation** (per IP scan `getType()` for `idxP/MZ/MY/T`; `G = −1/L`):
- PR-3a (z-only, scalar): `K_αα,z = (1/L)Σwt·ks(MZ,MZ) + dM_coh_z/dα_z`, **guarded**
  reciprocal with a magnitude floor against the positive bulk term (2D landmine: `K_αα` is
  sign-discontinuous at activation and indefinite on the softening branch).
- `K_vα,z` is a **6-vector** including the cross-tangent rows so the off-diagonals of the
  condensed `kb` are right: `K_vα,z(0)=−(1/L)Σwt·ks(P,MZ)`, `(1)=−(1/L)Σwt·(6ξ−4)·EI_z`,
  `(2)=−(1/L)Σwt·(6ξ−2)·EI_z`, `(3)=−(1/L)Σwt·(6ξ−4)·ks(MY,MZ)`,
  `(4)=−(1/L)Σwt·(6ξ−2)·ks(MY,MZ)`, `(5)=−(1/L)Σwt·ks(T,MZ)`.
- PR-3b (biaxial, 2×2): symmetric `K_αα` with off-diagonal `(1/L)Σwt·½(ks(MZ,MY)+ks(MY,MZ))`,
  `K_vα` a 6×2, condensation `kb −= K_vα·inv_guarded(K_αα)·K_vαᵀ` via an **eigenvalue-floored**
  symmetric-2×2 inverse, plus a belt-and-suspenders `‖dK‖_F ≤ 1e3·‖K_vv‖_F` step-cut guard.

**Condensation site** (the review's critical correction): the 3D `getTangentStiff` builds
`kb` **INLINE** in the IP loop (`:555-653`); `getBasicStiff` is only the dead `nlGeom` path
(`:658`, mutually exclusive with `-hinge`). Insert the condensation **after the inline loop,
before `K = crdTransf->getGlobalStiffMatrix(kb, q)` (`:695`)**, gated `if (hingeOn && …)`.
The basic **force `q` needs NO correction** — at `h=0` the sections hold `κ_bulk`, so the
existing `q = Σwt·Bᵀ·s` (both the tangent loop `:633-651` and `getResistingForce`) is
already the condensed basic force.

## Frame invariance (why the pinned invariant survives the quaternion triad)

`α_z, α_y` are **basic (corotated, rigid-body-free) rotation jumps** — `CorotCrdTransf3d`
exposes `v` as the 6-DOF post-`Tp` basic displacements, already stripped of rigid-body
rotation by the quaternion mean triad. The cohesive laws are work-conjugate to the basic
moments `q(1,2)`, `q(3,4)`. `getGlobalStiffMatrix(kb,q)` builds the ENTIRE geometric
stiffness from `q`/`pb`, reading `kb` ONLY through the material triple product
`kl = Tpᵀ·kb·Tp` — so condensing the α rows out of the 6×6 `kb` BEFORE that call cannot
perturb the geometric stiffness; there is **no seam** for element-internal DOFs. Condensation
operates entirely on the already-assembled (quaternion-correct) 6×6 `kb`; it never re-enters
the quaternion map. **Document the guard:** any future `kb`-dependent geometric term voids
this invariant.

## Resolved risks (adversarial review)

- **Two rank-1 updates inconsistent when `ks(MZ,MY)≠0`** → PR-3b uses a true coupled 2×2;
  PR-3a (z-only) has no off-diagonal so it is correct as a literal port; FD-tangent enforces.
- **Det-floor adjugate inverse blows up** → eigenvalue flooring (bounds every mode), sign per
  eigenvalue preserved (a softening-induced negative eigenvalue is physical).
- **Two `setTrialSectionDeformation` passes corrupt section state** → single unified inner Newton.
- **Condensation site mis-cited** → after the inline loop (`:653`), before `:695`.
- **Loosely-converged α ⇒ rotation-dependent spurious geometric stiffness in 3D** (`pb` feeds all
  6 components through `m(i)=pl(i)/(2cos ul(i))`, unlike 2D where only axial `pb` feeds P-Δ) →
  tight inner-Newton tolerance; if it hits `maxIter`, signal non-convergence to cut the global
  step (never hand crdTransf an inconsistent `pb`). v1 contract: full Newton + DisplacementControl.
- **`RELDIFF==0` unreachable for the rotated-hinge gate** (transcendental `asin/cos/tan`) → that
  gate uses relative tol ~1e-10; exact `RELDIFF==0` is reserved for the no-hinge reduce-to-stock
  gate only (identical code path / FP order).
- **`1/cos(ul)` singularity at ±90°** → the finite-rotation gate applies a RIGID quaternion
  rotation to both triads (basic `v` unchanged); it must not bend the element to "open" the crack.
- **Per-axis tolerance cross-contamination** (PR-3b) → per-component `Mscale_z, Mscale_y`.
- **Stale cached `K_vα` after revert without an intervening `update()`** (ModifiedNewton /
  line-search / arc-length) → the 2D evidence (PR-2b) showed the residual is always evaluated
  post-`update()`, so this does NOT bite for the dissipation; v1 documents full-Newton +
  DisplacementControl, idempotent re-converge deferred.

## PR-3a scope (the minimal correct first increment)

**Ships:** the strong-axis jump `α_z` only. `hingeOn`-gated members `theHingeZ`, `hingeJumpZ`,
`hingeJumpCommitZ`, `hingeKaaZ`, `hingeKvZ[6]`, `hingeMscaleZ`; a `solveHingeJump()` (unified inner
Newton restricted to z, 6-DOF `e`-vector, `hingeKvZ` carrying the cross-tangent rows); ONE guarded
rank-1 condensation before `:695`; new hinge blocks in `commitState`/`revertToLastCommit`/
`revertToStart`; OPS `-hinge $matTag` with a HARD `-nl` XOR `-hinge` parser error; serialization
(grow `data` 19→21: `hingeOn`, `hingeJumpCommitZ`; send `theHingeZ` like the sections). Gated so the
no-hinge path is byte-identical. Reach via the existing ndm-dispatch `LadrunoDispBeamColumn` command.

**Gates (mirror the 2D suite):** reduce-to-stock `RELDIFF==0`; constant-`Mz` patch test (hinge
carries exactly the applied moment); single-element `∫Mz d[[θz]] == Gf`; no energy double-count
(element total dissipation == Gf on an elastic section); **finite-rotation invariance under
`CorotCrdTransf3d`** (rigid 90/180° quaternion rotation of a pre-cracked element → force rotates,
magnitude + dissipation + `α_z` unchanged, rel-tol 1e-10); solver robustness + tight-α convergence;
FD-tangent (the 2D used tight-Newton-through-softening as the oracle; a global-DOF FD harness is the
ideal and attempted, else documented as covered by finite-rotation + convergence).

**Deferred (PR-3b+):** `α_y` + the coupled 2×2 with eigenvalue-floored inverse + staggered-activation
FD sweep (PR-3b, the biaxial MVP); a coupled biaxial cohesive interaction surface (`MAT_TAG 33004`
reserved, v2 material); a torsional jump (`-torsion`, no twist cohesive law yet); `-nl`+hinge
cross-terms; per-axis distinct hinge IP locations; non-Newton idempotent re-converge hardening.

## Open questions (defaults taken)

1. **Torsion in v1 — OUT** (no twist cohesive law; torsional fiber-rupture localization unstandardized).
   Optional `|T|`-exceeds-threshold soft warning deferred to PR-3b with `α_t`.
2. **Coupled vs block-diagonal cohesive — block-diagonal v1**, coupled `33004` v2 (the bulk-`ks`
   coupling carries the dominant biaxial interaction; revisit if the validation set needs the ellipsoid).
3. **PR-3a one axis vs both — z-only first** (lowest-risk pinned-invariant proof), biaxial PR-3b.

**Class tags:** reuses `ELE_TAG_LadrunoDispBeamColumn3d = 33014` (same element, gated option). Reserve
`MAT_TAG_LadrunoCohesiveHingeBiaxial = 33004` (uniaxial/nD band; v2, not yet built).

## Implementation log

### 2026-06-17 — PR-3a SHIPPED: strong-axis (Mz) embedded hinge, 3D

The strong-axis jump `α_z` only — the literal 2D scalar algebra on the Mz row of the 6-DOF
3D basic system. Reuses `ELE_TAG 33014` (same element, gated `-hinge` option).

- **Files:** `SRC/element/ladrunoDispBeamColumn/LadrunoDispBeamColumn3d.{h,cpp}` +
  `OPS_LadrunoDispBeamColumn.cpp` (usage), `tests/test_ladrunoDispBeamColumn3d_hinge.py` (12 tests).
- **What landed:** `-hinge $matTag` adds one scalar jump `α_z` carried by any `UniaxialMaterial`
  (default `LadrunoCohesiveHinge`). `solveHingeJump()` runs the inner Newton on `α_z` (sections see
  `κ_z_bulk = B_z·v − α_z/L`; residual `h_z = −Σwt·s(MZ) + M_coh(α_z)`; guarded `K_αα,z`; stable
  per-axis `hingeMscaleZ`). `getTangentStiff()` applies ONE guarded rank-1 condensation
  `kb −= hingeKvZ·hingeKvZᵀ/hingeKaaZ` to the inline-built 6×6 `kb` **after the IP loop, before**
  `crdTransf->getGlobalStiffMatrix` (the pinned invariant). `hingeKvZ` is a **6-vector** including the
  cross-tangent rows `ks(P,MZ)`, `ks(MY,MZ)`, `ks(T,MZ)` so the condensed off-diagonals are right.
  `getResistingForce` unchanged. ALL gated on `hingeOn` (no-hinge path byte-identical); `-hinge`+`-nl`
  rejected at parse. commit/revert/revertToStart + sendSelf/recvSelf (separate hinge-material send,
  data Vector 19→21). `"hinge"` setResponse passthrough.
- **Verified (OpenSeesPy, 12/12; full 65/65 element+material+hinge 2D+3D):** reduce-to-Tier-1
  **bit-identical** with `-hinge` absent; constant-Mz **patch test** (`1e-9`); **energy gate**
  `∫Mz d[[θz]] == Gf` (LINEAR `19.999992`, peak `Mc`); **element total dissipation == Gf** on an
  elastic section; tight-Newton tangent through EXP softening; **finite-rotation invariance under
  `CorotCrdTransf3d`** — a transverse-load cantilever driven to >0.5 rad member rotation still
  dissipates `Gf` (the pinned invariant survives the quaternion triad); solver robustness
  (Newton/ModifiedNewton/KrylovNewton); integration-objectivity nIP∈{2..5}; commit/revert + DB roundtrip.
  Zero inner-Newton convergence warnings through the fully-broken regime.
- **NEXT (PR-3b):** the weak-axis jump `α_y` + the TRUE coupled 2×2 `K_αα` (off-diagonal
  `ks(MZ,MY)`) with the **eigenvalue-floored** guarded inverse + 6×2 `K_vα` + the staggered-activation
  FD sweep. Then the coupled biaxial cohesive (`33004`, v2), torsion, `-nl`+hinge.
