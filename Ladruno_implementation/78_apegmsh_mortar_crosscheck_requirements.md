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
> (38 gates — mirror of the C++ kernel) and the C++ itself
> (`SRC/domain/contact/LadrunoMortarKernel.h`, `LadrunoTie.cpp::generateMortar`).
>
> **Revision 2 (2026-08-04)** — closes the loop with apeGmsh PR #898 (ADR 0086 S1
> merged; R7 case A reproduced to 5.0e-13). apeGmsh's report found four defects in
> THIS document, all fixed here: R7's master node numbering was misdescribed
> (§R7), R7's reference row was quoted too coarsely to support its own tolerance
> (§R7), R3 had no teeth because every R7 case was affine (§R3 + the new R7 case
> B), and R5.3's gap sub-gate is dead code in a flat-and-coincident scope (§R5.3).
> R8's ambiguity is resolved in §R8. R2 is unchanged and stays as written.

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
* Patch exactness is quadrature-independent (D and M share Gauss points; oracle
  T9 and T12 confirm it on affine AND non-affine facets). It is **P itself**,
  the guards' area/gap measures, and D conditioning that move with the rule.
  What "match the rule" actually buys depends on the facets (oracle T13):
  * **Insufficient degree is caught everywhere, affine included.** The degree-4
    rule under-integrates the degree-5/6 serendipity `N·φ` products even on
    rectangles: on R7 case A it shifts D by 6.4e-06 and `P_dual` by 4.7e-04.
    Use the degree-6 rule for quadratic pairs — that part of R3 has always had
    teeth.
  * **Two rules that are BOTH exact to degree 6 are indistinguishable on affine
    facets.** On rectangles and triangles the pullback integrand *is* a
    polynomial of degree ≤ 6, so Dunavant-12 and apeGmsh's Duffy-collapsed 5×5
    tensor rule (25 pts) agree to **1.9e-15** in `P_dual` — measured fork-side
    with apeGmsh's rule mirrored into the oracle (apeGmsh measured 4.5e-15).
    A porter can therefore ignore this whole section and still pass an
    affine-only R7. That was a real hole in revision 1.
  * **On non-affine facets neither rule is exact** — the isoparametric pullback
    of a non-parallelogram is rational, not polynomial — and the same two rules
    part company at **~1e-7** in `P_dual`, five orders above R7's tolerance.
    Fork-measured: 1.4e-07 on R7 case B, 3.5e-07 on a stronger skew; apeGmsh
    measured 8.4e-08 / 5.9e-08 on its own two skews. Independent, same verdict.
  ⇒ **Match the rule, and cross-check on R7 case B**, which is the only case that
  can enforce it. If you consciously diverge on the rule (apeGmsh v1 does), then
  say so and treat R7 case B as a **≤1e-7** agreement check rather than ≤1e-12;
  case A stays at ≤1e-12 either way.

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
   **Scope note (apeGmsh PR #898):** the mean-gap sub-gate only ever fires for
   *curved / non-coincident* master surfaces. In a **flat-and-coincident scope**
   — apeGmsh v1, and any port with the same restriction — every master already
   lies in one plane, so any finite mutual clip is a genuine overlap and the gap
   test can never veto the refusal. Implement the area test; the gap branch is
   dead code there and may be omitted (say so in the port's docs, per the
   "consciously diverge and document" rule above). Ship it only when the port
   admits curved or offset master surfaces.
4. Zero total overlap area ⇒ named refusal.

## R6 — things that stay fork-side (do NOT port)

Mass checks (HRZ / negative-lumped-mass refusals), EQ-constraint emission, the
projection handler. But one **export-side requirement lands on apeGmsh**: when
apeGmsh emits hex20 elements whose faces are mortar-tied and the run is
explicit (`LadrunoProjection`), the element must be `LadrunoBrick20 -lumped`
(HRZ) — row-sum lumping gives negative corner masses, which the fork now
refuses by name.

## R7 — the Q4 cross-check protocol (the deliverable that closes ADR 0086 D2)

Run both kernels on the SAME patches and compare P **row-by-row**. Two cases:
**A (affine)** pins the algebra, **B (non-affine)** pins the quadrature — case A
alone provably cannot (R3, oracle T13), which is why revision 2 adds B.

**Do not key the comparison to node indices.** Read the reference from
`kinematic_tie_validation/mortar_crosscheck_reference.json` (written by every
oracle run), which carries `slave_coords` / `master_coords` / `P_dual` / `P_std`
for both cases at **round-trip-exact float precision**, and match rows and
columns by coordinate. The digits printed below are for eyeballing only.

* **Case A — affine.** slave = 2×2 quad8 mesh on [0,1]², z=0 (21 nodes,
  corners-first per R1); master = 3×3 quad4 mesh on [0,1]², z=0 (16 nodes);
  refDir = +z; degree-6 rule. Tolerance **≤ 1e-12**.
  * `‖P_dual‖_F = 3.9279786887725048`, `‖P_std‖_F = 4.0157267126004861`
  * `P_dual` row for the slave node at (0,0), keyed by **master (x, y)**:

    | master (x, y) | weight |
    |---|---|
    | (0, 0)     | ` 1.1193415637860067` |
    | (1/3, 0)   | `-0.09053497942386797` |
    | (1/3, 1/3) | `-0.11522633744855894` |
    | (0, 1/3)   | `-0.09053497942386766` |
    | (2/3, 0)   | ` 0.04526748971193395` |
    | (2/3, 1/3) | ` 0.05761316872427926` |
    | (1/3, 2/3) | ` 0.057613168724279656` |
    | (0, 2/3)   | ` 0.04526748971193422` |
    | (2/3, 2/3) | `-0.02880658436213997` |

    (Sanity check on any port: this row is a corner node of a symmetric patch,
    so it **must** be symmetric under x ↔ y — the (1/3,0) and (0,1/3) weights
    agree, as do (2/3,0)/(0,2/3) and (2/3,1/3)/(1/3,2/3). Revision 1 quoted this
    row against node indices and called the master numbering "row-major"; it is
    neither row-major nor safe to rely on — it is the order the oracle's
    `quad4_mesh` `nid()` closure first assigns each coordinate as it walks
    `for ey: for ex:`. Under row-major indexing the pairs above do not line up
    and the row reads as spuriously asymmetric.)
* **Case B — non-affine (the R3 enforcement case).** The same topology,
  bilinearly mapped onto the trapezoid with corners (0,0), (1,0), (0.8,1),
  (0.2,1): `x = u + 0.2v − 0.4uv`, `y = v`. Mesh lines stay straight and edge
  midpoints stay midpoints, so R1/R2 are untouched and the R2 midside guard
  still passes — but the elements stop being parallelograms, the pullback goes
  rational, and quadrature stops being exact. Tolerance **≤ 1e-12 with the R3
  rule; ≤ 1e-7 if the port knowingly uses a different degree-6-exact rule** (R3).
  * `‖P_dual‖_F = 3.903529029143594`, `‖P_std‖_F = 3.992957735606101`
  * `P_dual` row for the slave node at (0,0), keyed by master (x, y):
    `(0,0): 1.0332251569787347`, `(1/3,0): −0.07429326901096864`,
    `(0.3555…,1/3): −0.030658239376578623`, `(0.0666…,1/3): 0.02317207474178856`,
    `(2/3,0): 0.03505202942130657`, `(0.6444…,1/3): 0.019518329856644683`,
    `(0.3777…,2/3): 0.01192051423475555`, `(0.1333…,2/3): −0.009881734644127314`,
    `(0.6222…,2/3): −0.008054862201555473`
    (asymmetric, correctly — the trapezoid map is not symmetric under x ↔ y).
* Additional cross-check cases already gated fork-side (T10): quad4-slave on
  two tri6 masters covering [0,1]²; quad8-slave on a 3×3 quad8 master mesh;
  tri3-slave on quad8 masters — all linear-patch exact ≤ 3e-15 in both bases.
  Note these are affine too, so they gate the BASIS, not the rule.
* **Agreed numbers (recorded both sides, 2026-08-04 — this closes ADR 0086 D2
  for case A).** apeGmsh PR #898 vs the fork oracle, case A:

  | quantity | fork | apeGmsh | diff |
  |---|---|---|---|
  | `‖P_dual‖_F` | 3.927978688773 | 3.927978688773 | 5.0e-13 |
  | reference corner row | above | matches to all quoted digits | — |
  | row sum | 1.000000000000 | 1.000000000000 | — |

  Pinned apeGmsh-side as `test_q4_crosscheck_matches_fork_oracle`
  (`tests/test_mortar_kernel.py`); fork-side by the oracle's Q4 block plus T12
  (case B) and T13 (why case B exists). apeGmsh reached case A with a Duffy 5×5
  rule rather than R3's Dunavant-12 — legitimate on affine facets, and precisely
  the divergence case B now measures. **Case B is not yet cross-checked**; it is
  published here for the next port (and for apeGmsh, should it adopt R3's rule).
* apeGmsh's port also found two genuine defects in its own kernel against this
  contract — a missing R2 midside guard and a missing R5.3 master self-overlap
  refusal, the latter behaving exactly as the R5.3 MAJOR finding predicts (the
  masked slave strip is extrapolated, so Σw = 1 and the linear patch both stay
  clean). The contract earns its keep; keep porting against it.

## R8 — test-honesty requirements (so apeGmsh's suite can actually fail)

* Include a linear-patch assertion **through each master type** (quad4, tri6,
  quad8) — partition-of-unity alone cannot catch a permuted midside **basis**.
* Mutation-test the suite once: corrupt one serendipity midside function and
  confirm the patch tests fail; restore.

**Which permutation R8 means (asked by apeGmsh PR #898 — answer: the basis).**
Revision 1 said "a cyclic tri6-midside permutation keeps Σφ=1, passes every
internal gate, and leaves a 0.15 patch error" without saying *what* was permuted.
It means the **shape functions — the basis itself miswired in code**, i.e. slot
`k` of the tri6 φ array evaluating the function that belongs to slot `k+1`. The
object under test is the kernel's own basis wiring (ADR-78 adversarial MAJOR-3:
"the C++ tri6/quad8 MASTER basis was verified by nothing"), not the mesh data.
R2 is irrelevant to it, and a **mutation test is the only instrument** — a
geometry guard cannot see a code bug. Both readings measured fork-side against
the oracle's T10(a) tri6-master patch:

| reading | R2 midside guard | Σφ = 1 | R4 PoU / P rowsums | R5 coverage, gapL1 | linear patch |
|---|---|---|---|---|---|
| permuted **connectivity** (midside node tags reordered) | **refuses, status −2** | — | — | — | never runs |
| permuted **basis** (φ slots miswired) | passes (blind to it) | 2.2e-16 | 5.6e-17 / 1.3e-15 | 1.000000, 0.0 | **0.35 (dual), 0.38 (std)** |

So for the **connectivity** reading, revision 1's "passes every internal gate" is
**not accurate for the fork** — permuted midsides are no longer at their edge
midpoints, the R2 guard refuses the facet structurally (fork status −2; the
oracle's `assemble` raises), and the patch test never gets to run. That is
exactly apeGmsh's observed behaviour, and it is correct behaviour. Only the
**basis** reading survives every gate, which is the reading R8 was written for.
(Revision 1's "0.15" came from the review's own probe geometry; the magnitude is
what matters — an O(0.1) patch error sailing past every internal gate. On the
oracle's published T10(a) patch it measures 0.35/0.38, reproducible by swapping
`shape_full`'s tri6 midside block for a cyclic permutation.)

apeGmsh implemented both readings — connectivity permutation asserted refused,
plus a live basis-mutation test — which is a strict superset of R8 and the right
call. Keep both.
