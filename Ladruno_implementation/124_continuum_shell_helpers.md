# WP-124 — Shared Element-contract helpers for the continuum element shells

Revision 1 (scoping complete). Refactor candidate 2 of WP-120 (#855) R3. The inventory was built by a
read-only Opus subagent; the four findings marked **verified** were re-checked by running the command.
Everything else carries the subagent's evidence (file:line) and is marked as such.

Status: **scoping done; draft PR #859.** No production code changed. Implementation waits for WP-123's local
builds (one build machine) and for the owner's choice on the gap fixes (C below).

Scoped 2026-09-25. Branch `wp/124-continuum-shell-helpers`, cut from `ladruno` @ `bc5c33453`.

## Problem

WP-120 R3 family 1: eight fork continuum elements (LadrunoQuad, LadrunoCST, LadrunoLST, LadrunoCSTPair,
LadrunoBrick, LadrunoBrick20, BezierTri6, BezierTet10) share an `Element`-contract shell (1,002 exact /
2,365 renamed-identifier shared lines). Eight PRs had to fix the same defect in several copies:

| PR | Defect | Copies |
|---|---|---|
| #562 | Rayleigh P-clobber in `getResistingForceIncInertia` | Quad, CST, LST, CSTPair |
| #228 | `getInitialStiff` not returning the cached `*Ki` | Quad, CST |
| #670 | mass cache not invalidated in `recvSelf` | 8 (incl. SolidShell) |
| #683 | EAS inner-Newton warning spam | Brick, Quad |
| #588 | EAS degeneracy guard blind to axis collapse | Brick, Quad |
| #709 | unsymmetric plastic tangent symmetrised | BezierTet10, BezierTri6 |
| #224 | `setParameter` not forwarded to GP materials | BezierTet10, BezierTri6 |
| #852 | ground-motion inertia sign | BezierTet10, BezierTri6 |

## Inventory (2026-09-25)

Six shared parts × eight elements. `=X` means identical to X. The variant legend is in the inventory's
source notes (file:line per cell).

| Part | Quad | CST | LST | CSTPair | Brick | Brick20 | Tri6 | Tet10 |
|---|---|---|---|---|---|---|---|---|
| 1 `getResistingForceIncInertia` | V1a | V1b | =Quad | V1b′ | V3a | V3b | V2 | =Tri6 |
| 2 `addInertiaLoadToUnbalance` | A1 | A2 | =Quad | A2′ | A4a | A4b | A3 | =Tri6 |
| 3 `getInitialStiff` cache | K1 | =Quad | =Quad | =Quad | **K2** | K1 + recv-invalidate | K1 | K1 + degeneracy guard |
| 4 mass cache | shared `LadrunoMassCache` | none | =Quad | none | ad hoc `Mi` | ad hoc `M0` + ρ key | =Quad | =Quad |
| 5 parameters | P1 | =Quad | =Quad | **absent** | P3 | P3′ | P2 | =Tri6 |
| 6 responses | R1 | =Quad | =Quad | R1′ | R2a | R2b | R3 | =Tri6 |

Operation order of part 1, which decides bit-identity:
- **V1 (plane):** ((f − Q) + diag(M)∘a) + R.
- **V2 (Bezier):** ((f − Q) + M·a) + R. For a diagonal M this is bit-equal to V1, apart from signed zeros and
  non-finite accelerations: `Vector::addMatrixVector` adds the zero off-diagonals.
- **V3 (Brick):** ((f + Σ_gp Nρ N a) + R) − Q. Here the load is subtracted after Rayleigh, so V3 cannot join
  V1/V2.

Constraint: `Element::getRayleighDampingForces()` and the factors are **protected** (`Element.h:206-209`), so a
free helper cannot do the Rayleigh tail. The snapshot + Rayleigh step stays in each element.

## Gaps: fixes that never reached every copy

| # | Gap | Evidence | Status |
|---|---|---|---|
| C1 | **#228 never reached LadrunoBrick**: its std/bbar/finite `getInitialStiff` builds `Ki`, then returns the class-static `stiff` | `LadrunoBrick.cpp:731-732`, `LadrunoBrick.h:250` (`static Matrix stiff;`) | **verified** |
| C2 | **Bezier never chains to `Element::setResponse`/`getResponse`**: globalForce, dampingForce, dynamicForce and inertialForce return nothing; no `stiffInitial`; no canonical aliases | `grep -c "Element::setResponse\|Element::getResponse"` = 0 in BezierTri6.cpp and BezierTet10.cpp (Quad: 3) | **verified** |
| C3 | **CSTPair has no `setParameter`/`updateParameter`**: `parameter`/`addToParameter` are a silent −1 (the class of defect #224 fixed for Bezier) | grep count 0 (CST: 4) | **verified** |
| C4 | CSTPair has no `stressPlaneStrain` (ID 21) token; the other three plane elements do | grep count: CST/LST/Quad 1, CSTPair 0 | **verified** |
| C5 | Quad SSP: `setParameter "material" k` targets slot k−1, but under SSP only slot 0 is live (setResponse and Brick remap) | Quad:1667 vs 1557, Brick:4243 | suspected (subagent) |
| C6 | `Ki` is invalidated in `recvSelf` only by Brick20, while Quad/LST/Brick/Tri6/Tet10 rewrite Ki-relevant state there | Brick20:1467; Quad:1422, LST:758, Brick:3611-3620, Tri6:1264, Tet10:1662-1673 | suspected: live-object re-receive only |
| C7 | Brick `updateParameter` keeps the last material's return; Brick20 fixed that ("F7") | Brick:4266-4270 vs Brick20:1305-1313 | divergence verified by subagent; reachability unverified |
| C8 | Brick `massType==1`: consistent residual inertia vs lumped `getMass` (tangent, αM, ground load). Brick20 fixed this pattern (F-1); Brick documents it as intentional | Brick:915 vs 926-929, 533-538; Brick20:862-912 | suspected |
| C9 | Brick's ad hoc mass-cache guard dereferences node pointers without the null check the shared helper gained | Brick:551 vs `LadrunoMassCache.h:87-90` | low impact |
| C10 | `materialState` token not excluded from the `material` branch in Quad/CST/LST/Brick/Brick20 (Bezier, SixNodeTri and LadrunoUP exclude it) | Tri6:1843, `SixNodeTri.cpp:1258`, LadrunoUP:1946 | suspected |
| C11 | Tri6 `getInitialStiff` caches even on a degenerate Jacobian; Tet10 refuses | Tri6:555-598 vs Tet10:476-480 | suspected |
| C12 | Bezier `recvSelf` ignores a material class change | Tri6:1286, Tet10:1691 | latent |
| C13 | Brick/Brick20: no `getRV` size check; Bezier reads `argv[0]` before an argc guard | Brick:787-789, Brick20:614-616, Tri6:1454, Tet10:1806 | minor |

## Shape (proposed)

Pure refactors, each stage proven bit-identical (fingerprint before/after, as WP-123) with a test that fails
on a one-line mutation of the helper. **Gap fixes are separate, labelled behaviour commits**, never mixed into
a refactor commit.

1. **Ki cache helper** (`cacheKi(Matrix *&Ki, const Matrix &formed)` + `dropKi`): Quad, CST, LST, CSTPair,
   Brick20, Tri6 first (already correct: a pure refactor), then Tet10. Then **C1** (Brick) as a fix commit,
   and C6 (drop Ki in `recvSelf`) as a fix commit.
2. **Ground-inertia helper** (`addGroundInertia(Q, M, nodes, nen, ndf, diagOnly, accel, who)`): Quad/LST
   (diagonal), then CST/CSTPair once they consume `getMass()`'s return value (retires the bare side-effect
   idiom), then Bezier/Brick/Brick20 (full M).
3. **Parameter and response helpers** (material forwarding with the "last non −1" rule; finalise =
   `endTag()` + qualified `Element::setResponse`). Fix commits for C2, C3, C4, C10 (and C5 after checking it).
4. **Brick mass cache → `LadrunoMassCache`** (C9). Brick20 keeps its own ρ-keyed design.
5. **Inertia-residual helper**, plane four + Bezier only (V1/V2). Brick's V3 stays; only its identical tail
   (Brick:819-826 ≡ Brick20:641-648) can be shared.

Riskiest: part 1 (GRFII), because of the per-element massless predicate and Rayleigh gate (V1 skips alphaM on
the massless path; V2 does not), diagonal vs full M·a, and the snapshot ordering (LEDGER_quirks "MUST snapshot
the shared static `resid`"). Do it last.

## Rejected approaches

- **A common base class.** The eight derive from `Element` with different node counts, DOF layouts and
  mass/stiffness caches, so a base would force one layout or a template hierarchy. Also, the Rayleigh tail
  needs `Element`'s protected members, which a free helper cannot reach but a base could; that is not worth
  a hierarchy change. WP-123's base worked because those four elements share a trivial, layout-free contract.
- **Folding the gap fixes into the refactor.** Each gap changes results. Mixed in, it would break the
  bit-identity proof that makes the refactor safe to merge.
- **Unifying Brick's V3 residual with V1/V2.** It subtracts the load after Rayleigh and accumulates inertia
  per Gauss point: not bit-identical.

## Open questions (owner)

- Fix the four verified gaps (C1–C4) in this WP as separate commits, or as their own small WPs? C1 and C2 are
  user-visible. For example, `recorder Element … globalForce` on a Bezier element currently records nothing.
- C8 (Brick lumped/consistent mismatch): documented as intentional in Brick; does the owner accept Brick20's
  F-1 reasoning for Brick too?
