# ADR-78 → apeGmsh: requirements for the numpy mortar kernel (ADR 0086 S1) and the Q4 cross-check

> Companion to fork ADR-78 (`78_ladruno_quadratic_mortar_tie_adr.md`, merged via
> PR #694). Audience: the apeGmsh team implementing
> `src/apeGmsh/_kernel/resolvers/_mortar.py` per apeGmsh ADR 0086
> (`0086-mortar-tie-composed.md`, D2 resolved "both"). The fork kernel is the
> REFERENCE implementation; this document is the contract the apeGmsh port must
> satisfy so the two kernels can cross-check row-by-row, plus the traps the
> fork's adversarial gate found so apeGmsh does not rediscover them.
> Working reference implementations (pure numpy, no build needed):
> `Ladruno_implementation/kinematic_tie_validation/proto_p2_2_quad8_mortar.py`
> (31 gates — mirror of the C++ kernel) and the C++ itself
> (`SRC/domain/contact/LadrunoMortarKernel.h`, `LadrunoTie.cpp::generateMortar`).

## R1 — facet vocabulary and node ordering (MUST match exactly)

* Supported facets: tri3, quad4, tri6, quad8. **Serendipity, corners-first**:
  tri6 = corners 0-2 then midsides 3-5 (midside k on edge (k, (k+1)%3));
  quad8 = corners 0-3 then midsides 4-7 (midside k on edge (k, (k+1)%4)).
  This is the gmsh / OpenSees 20-node-brick face convention — apeGmsh's
  `boundary_faces_for` extraction must deliver exactly this ordering.
* **tri6 SLAVE facets are a named refusal in BOTH bases** (corner ∫N = 0 over
  the facet ⇒ the dual scaling divides by zero — structural, no tolerance
  rescues it; the rowsum machinery is equally degenerate). tri6 MASTERS are
  fully supported. quad8 is fine both sides (corner ∫N = −A/12 ≠ 0).
* Basis functions (evaluate the FULL basis at the mapped Gauss point):
  quad8 corners ¼(1+ξξᵢ)(1+ηηᵢ)(ξξᵢ+ηηᵢ−1), mids ½(1−ξ²)(1±η) / ½(1∓ξ)(1−η²)
  per the ordering above; tri6 corners λ(2λ−1), mids 4λₐλᵦ with
  λ = (1−ξ−η, ξ, η).

## R2 — geometry: corner polygon only (the D1 design)

* ALL geometry runs on the **corner sub-facet** (first `ncorner = nps≤4 ? nps :
  nps/2` nodes): auxiliary plane from the master corners (centroid + corner
  normal, quad normal = cross of diagonals), Sutherland–Hodgman clip of the two
  corner polygons on that plane, fan triangulation from the clip centroid,
  inverse isomap of each Gauss point to slave corner parametrics, closest-point
  projection onto the master corner facet. Only the BASIS evaluation sees the
  full node set.
* Why exact: with midsides at edge midpoints the serendipity map reduces
  identically to the corner map — **including for warped (non-coplanar) corner
  quads** (review-verified to ~1e-15). Therefore:
  * **Guard (MUST): midside-at-midpoint** — `|X_mid − (X_a+X_b)/2| ≤
    1e-6 · max corner edge` (3D distance, catches curved edges). A violation is
    a **hard error naming the facet**, never a silent skip (the fork kernel
    returns a distinct status −2 vs −1 = "empty overlap, skip").
  * **No corner-coplanarity guard** — a warped quad8 degrades exactly like the
    (accepted) warped quad4: ~0.7 % constant-Jacobian area bias, patch tests
    unaffected. Do not add one; do not claim one.
* Projection Newton MUST carry the **scale-free step-converged escape**
  (`|dξ|+|dη| < 1e-8` after the update) alongside the absolute residual test —
  without it, models with coordinates ≳2000 length units silently lose most
  Gauss points (the fork's 2026-07 regression; its first oracle draft
  reproduced the bug by omitting it).

## R3 — quadrature

* Linear-linear pairs: the 6-point degree-4 rule (fork constants). Pairs with
  ANY quadratic facet: the **12-point degree-6 Dunavant rule** (serendipity
  N·φ products reach degree 6). Constants (barycentric, weights sum to 1):
  * w=0.116786275726379 at perms of (0.501426509658179, 0.249286745170910, ·);
  * w=0.050844906370207 at perms of (0.873821971016996, 0.063089014491502, ·);
  * w=0.082851075618374 at the 6 perms of (0.053145049844817, 0.310352451033784, ·).
  Self-verify exactness on monomials p+q ≤ 6 (oracle T2 does).
* Patch exactness is quadrature-independent (D and M share Gauss points), but
  the guards' area/gap measures and D conditioning are not — match the rule so
  the cross-check matches to 1e-12, not 1e-6.

## R4 — condensation

* Standard basis: global `P = D⁻¹M`. Dual (the production path; apeGmsh does
  not expose a toggle per ADR 0086): per-slave-facet
  `cᵉ_a = Σ_b Dᵉ_ab`, `Mdual += diag(cᵉ)·(Dᵉ)⁻¹Mᵉ`, `Ddual += cᵉ`,
  `P = Mdual / Ddual`.
* **quad8 corner dual masses are NEGATIVE (−A/12 tributary) by construction.**
  The sign cancels in P; any `Ddual > 0` assumption is a bug. Refusal must be
  sign-aware AND **relative**: `|Ddual[I]| ≤ 1e-12 · max_J|Ddual[J]|` (an
  absolute-zero test passes near-cancelled masses whose huge P rows the
  partition-of-unity check provably cannot catch — dual rowsums are identically
  1).
* Post-solve backstop (both bases): every P row sums to 1 within 1e-6, else
  refuse.

## R5 — guards (sign-free; the adversarial-gate lessons)

Per-node rowsum guards are meaningless for serendipity (negative corner
rowsums). The fork ships, and apeGmsh should mirror (or consciously diverge
and document):

1. **Coverage / protrusion, per slave facet**: Σ pair-clip areas vs the
   self-clip full facet area; refuse below `1 − 1e-3`.
2. **Conforming gap, per slave facet**: the kernel accumulates
   `gapL1 = ∫|g_N| dΓ` per pair; refuse `gapL1/areaCov > tol · 0.5·sqrt(areaFull)`
   (the 0.5 keeps the shipped per-node tributary strictness — a bare sqrt(A) is
   2× looser; found in review). The L1 form has no cancellation blind spot: an
   accordion surface whose shared-node SIGNED gaps cancel to ~1e-20 reads d/2
   in L1 (fork oracle T11).
3. **Master self-overlap refusal (MAJOR finding — do not skip)**: the coverage
   sum counts MULTIPLICITY, so a doubly-listed/overlapping master facet can
   exactly mask an uncovered slave strip. Before assembly, refuse any master
   pair with finite mutual clip area (> 1e-6 · facet area) whose mean gap is
   within the conforming tolerance. The gap gate keeps curved master surfaces
   legal (their aux-plane shadow overlaps have large gaps); edge-adjacent
   facets clip to slivers and drop out naturally.
4. Zero total overlap area ⇒ named refusal.

## R6 — things that stay fork-side (do NOT port)

Mass checks (HRZ / negative-lumped-mass refusals), EQ-constraint emission, the
projection handler. But one **export-side requirement lands on apeGmsh**: when
apeGmsh emits hex20 elements whose faces are mortar-tied and the run is
explicit (`LadrunoProjection`), the element must be `LadrunoBrick20 -lumped`
(HRZ) — row-sum lumping gives negative corner masses, which the fork now
refuses by name.

## R7 — the Q4 cross-check protocol (the deliverable that closes ADR 0086 D2)

Run both kernels on the SAME patch and compare **P row-by-row to ≤ 1e-12**:

* Geometry: slave = 2×2 quad8 mesh on [0,1]², z=0 (21 nodes, corners-first per
  R1); master = 3×3 quad4 mesh on [0,1]², z=0 (16 nodes); refDir = +z;
  degree-6 rule.
* Reference numbers (fork oracle, regenerate with
  `python proto_p2_2_quad8_mortar.py`):
  * `‖P_dual‖_F = 3.927978688773`, `‖P_std‖_F = 4.015726712600`;
  * `P_dual` row for slave node 0 at (0,0):
    `{0: 1.119341564, 1: −0.090534979, 2: −0.115226337, 3: −0.090534979,
      4: 0.045267490, 5: 0.057613169, 8: 0.057613169, 9: 0.045267490,
      10: −0.028806584}` (master node numbering = the oracle's `quad4_mesh`
    row-major order).
* Additional cross-check cases already gated fork-side (T10): quad4-slave on
  two tri6 masters covering [0,1]²; quad8-slave on a 3×3 quad8 master mesh;
  tri3-slave on quad8 masters — all linear-patch exact ≤ 3e-15 in both bases.
* Record the agreed numbers in BOTH repos (fork: extend the oracle's Q4 block;
  apeGmsh: the S1 unit tests). Until then the fork oracle carries the named Q4
  TODO.

## R8 — test-honesty requirements (so apeGmsh's suite can actually fail)

* Include a linear-patch assertion **through each master type** (quad4, tri6,
  quad8) — partition-of-unity alone cannot catch a permuted midside basis
  (review-probed: cyclic tri6-midside permutation keeps Σφ=1, passes every
  internal gate, and leaves a 0.15 patch error).
* Mutation-test the suite once: corrupt one serendipity midside function and
  confirm the patch tests fail; restore.
