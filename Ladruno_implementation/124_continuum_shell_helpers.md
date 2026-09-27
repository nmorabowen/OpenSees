# WP-124 — Shared Element-contract helpers for the continuum element shells

Revision 2 (implemented). Refactor candidate 2 of WP-120 (#855) R3. The inventory was built by a
read-only Opus subagent; its findings were then verified by running them (see **Gaps** and **Results**).

Status: **stages 1–5 done, every refactor bit-identical; gap fixes C1–C6, C9, C10, C14 landed as separate
behaviour commits; C7, C8, C11–C13 remain open (owner decisions, see Gaps).** Draft PR #859.

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
| C1 | **#228 never reached LadrunoBrick**: its std/bbar/finite `getInitialStiff` builds `Ki`, then returns the class-static `stiff` | `LadrunoBrick.cpp:731-732`, `LadrunoBrick.h:250` (`static Matrix stiff;`) | ✅ **fixed** `6e8bc9033` — invisible in a default build (every consumer copies at once); it was a hole in the ADR-87 CONTINUUM gate (rows C1a/C1b) |
| C2 | **Bezier never chains to `Element::setResponse`/`getResponse`**: globalForce, dampingForce, dynamicForce and inertialForce return nothing; no `stiffInitial`; no canonical aliases | `grep -c "Element::setResponse\|Element::getResponse"` = 0 in BezierTri6.cpp and BezierTet10.cpp (Quad: 3) | ✅ **fixed** `db1eb6875` — `test_base_response_vocabulary[Bezier*]` ("globalForce records nothing" before) |
| C3 | **CSTPair has no `setParameter`/`updateParameter`**: `parameter`/`addToParameter` are a silent −1 (the class of defect #224 fixed for Bezier) | grep count 0 (CST: 4) | ✅ **fixed** `56fa1ab50` — `parameter … E`: displacement ratio 1.0 before, 2.0 after |
| C4 | CSTPair has no `stressPlaneStrain` (ID 21) token; the other three plane elements do | grep count: CST/LST/Quad 1, CSTPair 0 | ✅ **fixed** `e10b8d560` — length 0 before, 8 after |
| C5 | Quad SSP: `setParameter "material" k` targets slot k−1, but under SSP only slot 0 is live (setResponse and Brick remap) | Quad:1667 vs 1557, Brick:4243 | ✅ **verified + fixed** `6fbfbdc25` — `material 2\|3\|4 E`: ratio exactly 1.0 vs 1.875 for `material 1 E` |
| C6 | `Ki` is invalidated in `recvSelf` only by Brick20, while Quad/LST/Brick/Tri6/Tet10 rewrite Ki-relevant state there | Brick20:1467; Quad:1422, LST:758, Brick:3611-3620, Tri6:1264, Tet10:1662-1673 | ✅ **verified + fixed** `9fb453372` — checkpoint restore into the live domain: `-initial` never converged (factor 1 − K/K0 = −1) on six elements |
| C7 | Brick `updateParameter` keeps the last material's return; Brick20 fixed that ("F7") | Brick:4266-4270 vs Brick20:1305-1313 | open — element-level `updateParameter` is only reached for element-registered parameters, none here; left for the owner |
| C8 | Brick `massType==1`: consistent residual inertia vs lumped `getMass` (tangent, αM, ground load). Brick20 fixed this pattern (F-1); Brick documents it as intentional | Brick:915 vs 926-929, 533-538; Brick20:862-912 | **observed**, open (owner decision): an elastic `-lumped` transient does not converge under Newton to 1e-10 (linear convergence) |
| C9 | Brick's ad hoc mass-cache guard dereferences node pointers without the null check the shared helper gained | Brick:551 vs `LadrunoMassCache.h:87-90` | ✅ **fixed** `33030b203` — latent (no script path found); made the stage-4 swap bit-identical on every path |
| C10 | `materialState` token not excluded from the `material` branch in Quad/CST/LST/Brick/Brick20 (Bezier, SixNodeTri and LadrunoUP exclude it) | Tri6:1843, `SixNodeTri.cpp:1258`, LadrunoUP:1946 | ✅ **verified + fixed** `05524f3f8` — `materialState` unclaimed on the five (DruckerPrager) before; Bezier as control |
| C11 | Tri6 `getInitialStiff` caches even on a degenerate Jacobian; Tet10 refuses | Tri6:555-598 vs Tet10:476-480 | open |
| C12 | Bezier `recvSelf` ignores a material class change | Tri6:1286, Tet10:1691 | open |
| C13 | Brick/Brick20: no `getRV` size check; Bezier reads `argv[0]` before an argc guard | Brick:787-789, Brick20:614-616, Tri6:1454, Tet10:1806 | open — the ground-inertia helper keeps the bricks' no-check behaviour (`checkSize=false`) |
| C14 | **Brick20 singular after a live restore**: `recvSelf` cleared `geomCached` and relied on the `setDomain` that follows a broker-built receive; the live branch of `Domain::recvSelf` calls `recvSelf` + `update()` only | found in WP-124 (C6 probe) | ✅ **fixed** `4849f0c2f` — `test_live_restore_is_usable` (singular before, Brick20 only) |

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

As built: all five stages as planned. Two things were NOT shared: the Brick/Brick20 V3 residual tail (kept in
the bricks, their residual is not M·a) and the massless predicates (they differ per element and stay local).
Stage 4 added `LadrunoMassCache::cached()` so the brick keeps returning its per-instance copy on a miss.

## Results (2026-09-27)

**Proof method.** Baseline = the unchanged code (`ladruno` merged at the WP start). `fingerprint.py --suite
shells` (models: `wp124_shells/shell_models.py`): 35 formulation variants × static DisplacementControl into the
plastic range (LadrunoJ2) with every setResponse token, ModifiedNewton -initial, algorithm Linear, eigen,
Newmark with each Rayleigh factor ALONE from a plastic preload (K, Kc, K0 differ), HHT, UniformExcitation,
CentralDifference, element parameters. 707 series, 1,884,982 values; two baseline runs identical. Builds are
done from the uncommitted tree and committed after the evidence (a named `build.bat OpenSeesPy` is incremental
then; the first pyd build of a tree compiles everything).

| Step | Fingerprint vs previous step |
|---|---|
| stage 1 Ki cache (1a, 1b, 1c) + C1 + C6 | 0 / 707 differ |
| stage 2 ground inertia (2a, 2b, 2c) | 0 / 707 |
| C14 | 0 / 707 |
| stage 3 parameters + responses | 0 / 707 |
| C2, C3, C4, C5, C10 (one build) | 170 differ, **all explained**: 163 only inside the NEW response blocks (Bezier globalForce/dampingForce/dynamicForce/inertialForce, CSTPair stressPlaneStrain; every other value byte-equal once those blocks are stripped) + 7 CSTPair parameter runs whose parameter now reaches the material (C3); 0 unexplained |
| C9 + stage 4 Brick → LadrunoMassCache | 0 / 707 (against the post-fix fingerprint) |
| stage 5 inertia residual | 0 / 707 |

**Mutation rows** (`wp124_shells/mutation_rows.py`: one-line edit of the helper, pyd rebuilt, sources restored;
test file `tests/test_ladruno_element_shell_helpers.py`):

| Row | Mutation | Result |
|---|---|---|
| K1 | `cacheKi` stores a zero matrix | caught (21) |
| K3 | `dropKi` keeps the stale Ki | caught (6) |
| K4 | `cacheKi` returns the scratch (the C1 shape) | EQUIVALENT in a default build: every consumer copies at once (predicted) |
| C1a / C1b | CONTINUUM=IDENT mutant, C1 fix reverted / kept | probe fails (the mutation never reaches the solver) / passes |
| G1 / G2 | ground inertia sign, full / diagonal branch | caught (8 / 8) |
| G3 | full M reduced to its diagonal | caught (4 → 8 after adding `-cMass` Bezier probes: Bezier defaults to LUMPED) |
| G4 | every DOF takes the x ground acceleration | caught (16) |
| P1 / P2 | forall asks only material 1 / every `material k` → slot 0 | caught (7 / 7; CST has one GP) |
| P3 | `finishResponse` drops the `Element::setResponse` fallback | caught (8) |
| P4 | `materialState` taken as a GP address again | caught (5) |
| P5 | `finishResponse` does not `endTag()` first | EQUIVALENT: `eleResponse` writes no XML, and `XmlFileStream` closes open tags itself (an XML-recorder well-formedness test also passes); the rule stays documented, no observer found |
| M1 | Brick cache signature without ρ | caught (1) — only after the rho-guard test moved from `inertialForce` to `dampingForce`: Brick's residual integrates ρ per GP and never reads the cache, so the first version of the test could not see M1 |
| N1 / N2 | residual inertia sign, diagonal / full | caught (14 / 14) |
| N3 | full M reduced to its diagonal in the residual | caught (6) |
| N4 | every DOF takes the x trial acceleration | caught (24) |
| C14 | Brick20 `recvSelf` without the rebuild | caught (2) |

**Gap evidence** is in the Gaps table (each fix's commit and before/after).

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

- ~~Fix C1–C4 here?~~ Decided: here, as separate behaviour commits (done, with C5, C6, C9, C10, C14).
- C8 (Brick lumped/consistent mismatch): documented as intentional in Brick; does the owner accept Brick20's
  F-1 reasoning for Brick too? Now observed: an elastic `-lumped` transient converges only linearly under
  Newton (the fingerprint's `Brick/lumped` needs a 1e-7 / 300-iteration test to run at all).
- C7, C11, C12, C13 stay open (low impact; C13's brick size check is behind `checkSize=false`).
